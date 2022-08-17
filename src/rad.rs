use crate::annotation;
use crate::convert;

extern crate fnv;

use libradicl::rad_types::decode_int_type_tag;
//use libradicl::rad_types::write_str_bin;
//use libradicl::rad_types::encode_type_tag;
use std::error::Error;
use std::fs::{self, File};
use std::path::Path;

use crate::query_bam_records::BAMQueryRecordReader;
use annotation::ExonNode;

use coitrees::COITree;
use fnv::FnvHashMap;
use rust_htslib::{bam, bam::record, bam::record::Aux, bam::Read};
use std::io::{BufWriter, Cursor, Seek, SeekFrom, Write};
// use scroll::Pread;

// fn read_into_u64<T: std::io::Read>(reader: &mut T, rt: &libradicl::rad_types::RadIntId) -> u64 {
//     let mut rbuf = [0u8; 8];
//     let v: u64;
//     match rt {
//         libradicl::rad_types::RadIntId::U8 => {
//             reader.read_exact(&mut rbuf[0..1]).unwrap();
//             v = rbuf.pread::<u8>(0).unwrap() as u64;
//         }
//         libradicl::rad_types::RadIntId::U16 => {
//             reader.read_exact(&mut rbuf[0..2]).unwrap();
//             v = rbuf.pread::<u16>(0).unwrap() as u64;
//         }
//         libradicl::rad_types::RadIntId::U32 => {
//             reader.read_exact(&mut rbuf[0..4]).unwrap();
//             v = rbuf.pread::<u32>(0).unwrap() as u64;
//         }
//         libradicl::rad_types::RadIntId::U64 => {
//             reader.read_exact(&mut rbuf[0..8]).unwrap();
//             v = rbuf.pread::<u64>(0).unwrap();
//         }
//     }
//     v
// }

// from https://github.com/COMBINE-lab/alevin-fry/blob/master/libradicl/src/convert.rs#L63-L85
pub fn cb_string_to_u64(cb_str: &[u8]) -> Result<u64, Box<dyn Error>> {
    let mut cb_id: u64 = 0;
    for (idx, nt) in cb_str.iter().rev().enumerate() {
        let offset = idx * 2;
        match nt {
            65 => (),                   // A 00
            67 => cb_id |= 1 << offset, // C 01
            71 => cb_id |= 2 << offset, // G 10
            84 => cb_id |= 3 << offset, // T 11
            _ => panic!("unknown nucleotide {}", nt),
        };
    }

    Ok(cb_id)
}

fn bam_peek(bam_path: &String) -> record::Record {
    let mut bam_reader = bam::Reader::from_path(&bam_path).unwrap();
    bam_reader.set_threads(1).unwrap();
    //
    let mut rec = record::Record::new();
    let first_record_exists = bam_reader.read(&mut rec).is_some();
    if !first_record_exists {
        log::error!("bam file had no records!");
        std::process::exit(1);
    }
    rec
}

pub fn bam2rad_bulk(
    input_bam_filename: &String,
    output_rad_filename: &String,
    transcripts: &Vec<String>,
    txp_lengths: &Vec<i32>,
    trees: &FnvHashMap<String, COITree<ExonNode, u32>>,
    threads_count: &usize,
    max_softlen: &usize,
) {
    let first_record = bam_peek(input_bam_filename);
    ////////////////////////////////////////////////// header section
    // is paired-end
    if first_record.is_paired() {
        bam2rad_bulk_pe(
            input_bam_filename,
            output_rad_filename,
            transcripts,
            txp_lengths,
            trees,
            threads_count,
            max_softlen,
        );
    } else {
        bam2rad_bulk_se(
            input_bam_filename,
            output_rad_filename,
            transcripts,
            txp_lengths,
            trees,
            threads_count,
            max_softlen,
        );
    }
}

fn dump_collected_alignments_bulk_se(all_read_records: &Vec<record::Record>, owriter: &mut Cursor<Vec<u8>>) -> bool {
    // add stored data to the current chunk
    let mut wrote_some: bool = false;
    let mut written_records: u32 = 0;
    let mut data = Cursor::new(vec![]);

    for txp_rec in all_read_records.iter() {
        if !txp_rec.is_unmapped() {
            // alignment reference ID
            // data.write_all(&(txp_rec.tid() as u32).to_le_bytes()).unwrap();
            // alignment orientation
            // data.write_all(&(txp_rec.is_reverse() as u8).to_le_bytes()).unwrap();
            // alignment position
            // data.write_all(&(txp_rec.pos() as u32).to_le_bytes()).unwrap();

            // array of alignment-specific tags
            // compressed_ori_refid
            let mut tid_comressed = txp_rec.tid() as u32;
            if !txp_rec.is_reverse() {
                tid_comressed |= 0x80000000 as u32;
            }
            data.write_all(&tid_comressed.to_le_bytes()).unwrap();
            // alnscore
            // NOTE: since AS is of type integer but RAD format doesn't have signed integer, we store it in a floating point number instead
            let mut alnscore: f32 = 0.0;
            match txp_rec.aux(b"AS") {
                Ok(value) => match value {
                    Aux::I8(v) => alnscore = v as f32,
                    Aux::U8(v) => alnscore = v as f32,
                    Aux::I16(v) => alnscore = v as f32,
                    Aux::U16(v) => alnscore = v as f32,
                    Aux::I32(v) => alnscore = v as f32,
                    Aux::U32(v) => alnscore = v as f32,
                    _ => println!("something else"),
                },
                Err(e) => {
                    log::error!(
                        "Could not find AS tag for a record of the read {}",
                        String::from_utf8(txp_rec.qname().to_vec()).unwrap()
                    );
                    panic!("Error reading AS tag: {}", e);
                }
            }
            data.write_all(&alnscore.to_le_bytes()).unwrap();
            // alnpos
            let alnpos: u32 = txp_rec.pos() as u32;
            data.write_all(&alnpos.to_le_bytes()).unwrap();

            written_records += 1;
            wrote_some = true;
        }
    }

    if written_records > 0 {
        // number of alignments
        owriter.write_all(&written_records.to_le_bytes()).unwrap();
        // read-level tags
        owriter
            .write_all(&(all_read_records.first().unwrap().seq_len() as u16).to_le_bytes())
            .unwrap(); // read length
                       // dump records
        owriter.write_all(data.get_ref()).unwrap();
    }

    return wrote_some;
}

pub fn bam2rad_bulk_se(
    input_bam_filename: &String,
    output_rad_filename: &String,
    transcripts: &Vec<String>,
    txp_lengths: &Vec<i32>,
    trees: &FnvHashMap<String, COITree<ExonNode, u32>>,
    threads_count: &usize,
    max_softlen: &usize,
) {
    let ofile = File::create(&output_rad_filename).unwrap();
    // file writer and intermediate buffer
    let mut owriter = BufWriter::with_capacity(1048576, ofile);
    let mut data = Cursor::new(vec![]);
    //
    let is_paired: u8;
    let end_header_pos: u64;

    {
        ////////////////////////////////////////////////// header section
        // is paired-end
        is_paired = 0;

        data.write_all(&is_paired.to_le_bytes()).expect("couldn't write to output file");
        // number of references
        let ref_count = transcripts.len() as u64;
        log::info!("Number of reference sequences: {}", ref_count);
        data.write_all(&ref_count.to_le_bytes()).expect("couldn't write to output file");
        // list of reference names
        for tname in transcripts.iter() {
            let name_len = tname.len() as u16;
            data.write_all(&name_len.to_le_bytes()).expect("coudn't write to output file");
            data.write_all(tname.as_bytes()).expect("coudn't write to output file");
        }
        // keep a pointer to header pos
        // end_header_pos = data.seek(SeekFrom::Current(0)).unwrap() - std::mem::size_of::<u64>() as u64;
        end_header_pos = data.seek(SeekFrom::Current(0)).unwrap();
        log::debug!("end header pos: {:?}", end_header_pos);
        // number of chunks
        let initial_num_chunks = 0u64;
        data.write_all(&initial_num_chunks.to_le_bytes()).expect("coudn't write to output file");

        ///////////////////////////////////////// tag definition section
        // FILE-LEVEL tags
        let mut num_tags = 0u16; // no file-level tags for bulk
        data.write_all(&num_tags.to_le_bytes()).expect("coudn't write to output file");

        // READ-LEVEL tags
        num_tags = 1u16; // no read-level tags for bulk
        data.write_all(&num_tags.to_le_bytes()).expect("coudn't write to output file");

        // read length
        let mut tag_str = "readlen";
        let mut tag_typeid = 2u8;
        libradicl::rad_types::write_str_bin(&tag_str, &libradicl::rad_types::RadIntId::U16, &mut data);
        data.write_all(&tag_typeid.to_le_bytes()).expect("coudn't write to output file");

        // ALIGNMENT-LEVEL tags
        num_tags = 3u16;
        data.write_all(&num_tags.to_le_bytes()).expect("couldn't write to output file");

        // reference id
        tag_str = "compressed_ori_refid";
        tag_typeid = 3u8;
        libradicl::rad_types::write_str_bin(&tag_str, &libradicl::rad_types::RadIntId::U16, &mut data);
        data.write_all(&tag_typeid.to_le_bytes()).expect("coudn't write to output file");

        // alignment score
        tag_str = "alnscore";
        tag_typeid = 5u8;
        libradicl::rad_types::write_str_bin(&tag_str, &libradicl::rad_types::RadIntId::U16, &mut data);
        data.write_all(&tag_typeid.to_le_bytes()).expect("coudn't write to output file");

        // reference position
        tag_str = "alnpos";
        tag_typeid = 3u8;
        libradicl::rad_types::write_str_bin(&tag_str, &libradicl::rad_types::RadIntId::U16, &mut data);
        data.write_all(&tag_typeid.to_le_bytes()).expect("coudn't write to output file");

        ///////////////////////////////////////// file-level tag values
        // nothing for bulk!
    }

    // dump current buffer content
    owriter.write_all(data.get_ref()).unwrap();

    // let mut bam_reader = bam::Reader::from_path(&input_bam_filename).unwrap();
    // let header_view = bam_reader.header().to_owned();

    // if *threads_count > 1 {
    //     bam_reader.set_threads(threads_count - 1).unwrap();
    // } else {
    //     bam_reader.set_threads(1).unwrap();
    // }

    let reader_threads: Option<usize>;
    if *threads_count > 1 {
        log::info!("thread count: {}", threads_count);
        reader_threads = Some(threads_count - 1);
    } else {
        reader_threads = None;
        log::info!("thread count: {}", threads_count);
    }

    // setup the input BAM Reader
    let mut bqr = BAMQueryRecordReader::new(input_bam_filename, reader_threads);
    let input_header = bqr.get_header().to_owned();

    let mut total_num_chunks = 0u64;
    let mut chunk_reads = 0u32;

    // allocate buffer
    let buf_limit = 10000u32;
    data = Cursor::new(Vec::<u8>::with_capacity((buf_limit * 100) as usize));
    // placeholder for number of bytes and number of records
    data.write_all(&chunk_reads.to_le_bytes()).unwrap();
    data.write_all(&chunk_reads.to_le_bytes()).unwrap();

    let mut all_query_records: Vec<record::Record> = Vec::new();
    let required_tags: Vec<&str> = vec![];

    while let Ok(Some(ret_vec)) = bqr.get_next_query_records() {
        all_query_records.clear();
        for r in ret_vec.iter() {
            let mut txp_records = convert::convert_query_bam_records(r, &input_header, transcripts, txp_lengths, trees, max_softlen, &required_tags);
            all_query_records.append(&mut txp_records);
        }
        if !all_query_records.is_empty() {
            let write_success = dump_collected_alignments_bulk_se(&all_query_records, &mut data);
            if write_success {
                chunk_reads += 1;
            }
        }
        if chunk_reads >= buf_limit {
            // dump current chunk and start a new one
            data.set_position(0);
            // number of bytes
            data.write_all(&(data.get_ref().len() as u32).to_le_bytes()).unwrap();
            // number reads
            data.write_all(&chunk_reads.to_le_bytes()).unwrap();
            owriter.write_all(data.get_ref()).unwrap();
            total_num_chunks += 1;
            chunk_reads = 0;
            data = Cursor::new(Vec::<u8>::with_capacity((buf_limit * 100) as usize));
            // placeholder for number of bytes and number of records
            data.write_all(&chunk_reads.to_le_bytes()).unwrap();
            data.write_all(&chunk_reads.to_le_bytes()).unwrap();
        }
    }
    // dump the last chunk
    if chunk_reads > 0 {
        data.set_position(0);
        // number of bytes
        data.write_all(&(data.get_ref().len() as u32).to_le_bytes()).unwrap();
        // number reads
        data.write_all(&chunk_reads.to_le_bytes()).unwrap();
        owriter.write_all(data.get_ref()).unwrap();
        total_num_chunks += 1;
    }
    // update the number of chunks in the header
    owriter.flush().expect("File buffer could not be flushed");
    owriter.seek(SeekFrom::Start(end_header_pos)).expect("couldn't seek in output file");
    owriter
        .write_all(&total_num_chunks.to_le_bytes())
        .expect("couldn't write to output file.");
}

fn dump_collected_alignments_bulk_pe(all_read_records: &[record::Record], owriter: &mut Cursor<Vec<u8>>) -> bool {
    // add stored data to the current chunk
    let mut wrote_some: bool = false;
    let mut written_records: u32 = 0;
    let mut data = Cursor::new(vec![]);

    let mut pos_left: u32;
    let mut pos_right: u32;
    let mut len_left: u16 = 0;
    let mut len_right: u16 = 0;

    let mut rec_iter = all_read_records.iter();
    while let Some(rec1) = rec_iter.next() {
        // there should be another mate
        if let Some(rec2) = rec_iter.next() {
            if rec1.is_unmapped() && rec2.is_unmapped() {
                // both mates are unmapped; no need to store anything in RAD format
                continue;
            } else if !rec1.is_unmapped() && !rec2.is_unmapped() {
                // both mates are mapped
                // alignment reference ID
                assert_eq!(rec1.tid(), rec2.tid());
                data.write_all(&(rec1.tid() as u32).to_le_bytes()).unwrap();
                let mut aln_type: u8 = 0;
                if rec1.is_first_in_template() {
                    // rec1 is left
                    if rec1.is_reverse() {
                        aln_type |= 2u8;
                    }
                    if rec2.is_reverse() {
                        aln_type |= 1u8;
                    }
                    pos_left = rec1.pos() as u32;
                    pos_right = rec2.pos() as u32;
                    len_left = rec1.seq_len() as u16;
                    len_right = rec2.seq_len() as u16;
                } else {
                    // rec2 is left
                    if rec1.is_reverse() {
                        aln_type |= 1u8;
                    }
                    if rec2.is_reverse() {
                        aln_type |= 2u8;
                    }
                    pos_right = rec1.pos() as u32;
                    pos_left = rec2.pos() as u32;
                    len_right = rec1.seq_len() as u16;
                    len_left = rec2.seq_len() as u16;
                }
                // alignment type (0..7)
                data.write_all(&aln_type.to_le_bytes()).unwrap();

                // array of alignment-specific tags
                // alnscore
                // NOTE: since AS is of type integer but RAD format doesn't have signed integer, we store it in a floating point number instead
                let mut alnscore: f32 = 0.0;
                match rec1.aux(b"AS") {
                    Ok(value) => match value {
                        Aux::I8(v) => alnscore = v as f32,
                        Aux::U8(v) => alnscore = v as f32,
                        Aux::I16(v) => alnscore = v as f32,
                        Aux::U16(v) => alnscore = v as f32,
                        Aux::I32(v) => alnscore = v as f32,
                        Aux::U32(v) => alnscore = v as f32,
                        _ => println!("something else"),
                    },
                    Err(e) => {
                        log::error!(
                            "Could not find AS tag for a record of the read {}",
                            String::from_utf8(rec1.qname().to_vec()).unwrap()
                        );
                        panic!("Error reading AS tag: {}", e);
                    }
                }
                data.write_all(&alnscore.to_le_bytes()).unwrap();
                data.write_all(&pos_left.to_le_bytes()).unwrap();
                data.write_all(&pos_right.to_le_bytes()).unwrap();
                data.write_all(&(rec1.insert_size().unsigned_abs()).to_le_bytes()).unwrap();
                wrote_some = true;
                written_records += 1;
            } else {
                // only one mate is mapped
                let mapped_rec: &record::Record;
                let unmapped_rec: &record::Record;
                if rec1.is_unmapped() {
                    mapped_rec = rec2;
                    unmapped_rec = rec1;
                } else {
                    mapped_rec = rec1;
                    unmapped_rec = rec2;
                }
                // there is no mate
                // alignment reference ID
                data.write_all(&(mapped_rec.tid() as u32).to_le_bytes()).unwrap();
                let aln_type: u8;
                if mapped_rec.is_first_in_template() {
                    // right is unmapped
                    aln_type = if mapped_rec.is_reverse() { 5 } else { 4 };
                    pos_left = mapped_rec.pos() as u32;
                    pos_right = 0u32;
                    len_left = mapped_rec.seq_len() as u16;
                    len_right = unmapped_rec.seq_len() as u16;
                } else {
                    // left is unmapped
                    aln_type = if mapped_rec.is_reverse() { 7 } else { 6 };
                    pos_left = 0u32;
                    pos_right = mapped_rec.pos() as u32;
                    len_left = unmapped_rec.seq_len() as u16;
                    len_right = mapped_rec.seq_len() as u16;
                }
                // alignment type (0..7)
                data.write_all(&aln_type.to_le_bytes()).unwrap();

                // array of alignment-specific tags
                // alnscore
                // NOTE: since AS is of type integer but RAD format doesn't have signed integer, we store it in a floating point number instead
                let mut alnscore: f32 = 0.0;
                match mapped_rec.aux(b"AS") {
                    Ok(value) => match value {
                        Aux::I8(v) => alnscore = v as f32,
                        Aux::U8(v) => alnscore = v as f32,
                        Aux::I16(v) => alnscore = v as f32,
                        Aux::U16(v) => alnscore = v as f32,
                        Aux::I32(v) => alnscore = v as f32,
                        Aux::U32(v) => alnscore = v as f32,
                        _ => println!("something else"),
                    },
                    Err(e) => {
                        log::error!(
                            "Could not find AS tag for a record of the read {}",
                            String::from_utf8(mapped_rec.qname().to_vec()).unwrap()
                        );
                        panic!("Error reading AS tag: {}", e);
                    }
                }
                data.write_all(&alnscore.to_le_bytes()).unwrap();
                data.write_all(&pos_left.to_le_bytes()).unwrap();
                data.write_all(&pos_right.to_le_bytes()).unwrap();
                data.write_all(&(mapped_rec.insert_size().unsigned_abs()).to_le_bytes()).unwrap();
                wrote_some = true;
                written_records += 1;
            }
        } else {
            log::error!("couldn't find respective mate!");
            return false;
        }
    }

    if written_records > 0 {
        // number of alignments
        owriter.write_all(&written_records.to_le_bytes()).unwrap();
        // read-level tags
        // first mate length
        owriter.write_all(&len_left.to_le_bytes()).unwrap();
        // second mate length
        owriter.write_all(&len_right.to_le_bytes()).unwrap();
        // dump records
        owriter.write_all(data.get_ref()).unwrap();
    }

    wrote_some
}

pub fn bam2rad_bulk_pe(
    input_bam_filename: &String,
    output_rad_filename: &String,
    transcripts: &Vec<String>,
    txp_lengths: &Vec<i32>,
    trees: &FnvHashMap<String, COITree<ExonNode, u32>>,
    threads_count: &usize,
    max_softlen: &usize,
) {
    let ofile = File::create(&output_rad_filename).unwrap();
    // file writer and intermediate buffer
    let mut owriter = BufWriter::with_capacity(1048576, ofile);
    let mut data = Cursor::new(vec![]);
    //
    let is_paired: u8;
    let end_header_pos: u64;

    {
        ////////////////////////////////////////////////// header section
        // is paired-end
        is_paired = 1;

        data.write_all(&is_paired.to_le_bytes()).expect("couldn't write to output file");
        // number of references
        let ref_count = transcripts.len() as u64;
        log::info!("Number of reference sequences: {}", ref_count);
        data.write_all(&ref_count.to_le_bytes()).expect("couldn't write to output file");
        // list of reference names
        for tname in transcripts.iter() {
            let name_len = tname.len() as u16;
            data.write_all(&name_len.to_le_bytes()).expect("coudn't write to output file");
            data.write_all(tname.as_bytes()).expect("coudn't write to output file");
        }
        // keep a pointer to header pos
        // end_header_pos = data.seek(SeekFrom::Current(0)).unwrap() - std::mem::size_of::<u64>() as u64;
        end_header_pos = data.seek(SeekFrom::Current(0)).unwrap();
        log::debug!("end header pos: {:?}", end_header_pos);
        // number of chunks
        let initial_num_chunks = 0u64;
        data.write_all(&initial_num_chunks.to_le_bytes()).expect("coudn't write to output file");

        ///////////////////////////////////////// tag definition section
        // FILE-LEVEL tags
        let mut num_tags = 0u16; // no file-level tags for bulk
        data.write_all(&num_tags.to_le_bytes()).expect("coudn't write to output file");

        // READ-LEVEL tags
        num_tags = 2u16; // no read-level tags for bulk
        data.write_all(&num_tags.to_le_bytes()).expect("coudn't write to output file");

        // read length left mate
        let mut tag_str = "readlen_left";
        let mut tag_typeid = 2u8;
        libradicl::rad_types::write_str_bin(tag_str, &libradicl::rad_types::RadIntId::U16, &mut data);
        data.write_all(&tag_typeid.to_le_bytes()).expect("coudn't write to output file");

        // read length right mate
        tag_str = "readlen_right";
        tag_typeid = 2u8;
        libradicl::rad_types::write_str_bin(tag_str, &libradicl::rad_types::RadIntId::U16, &mut data);
        data.write_all(&tag_typeid.to_le_bytes()).expect("coudn't write to output file");

        // ALIGNMENT-LEVEL tags
        num_tags = 6u16;
        data.write_all(&num_tags.to_le_bytes()).expect("couldn't write to output file");

        // reference id
        tag_str = "refid";
        tag_typeid = 3u8;
        libradicl::rad_types::write_str_bin(tag_str, &libradicl::rad_types::RadIntId::U16, &mut data);
        data.write_all(&tag_typeid.to_le_bytes()).expect("coudn't write to output file");

        // alignment type
        tag_str = "alntype";
        tag_typeid = 1u8;
        libradicl::rad_types::write_str_bin(tag_str, &libradicl::rad_types::RadIntId::U16, &mut data);
        data.write_all(&tag_typeid.to_le_bytes()).expect("coudn't write to output file");

        // alignment score
        tag_str = "alnscore";
        tag_typeid = 5u8;
        libradicl::rad_types::write_str_bin(tag_str, &libradicl::rad_types::RadIntId::U16, &mut data);
        data.write_all(&tag_typeid.to_le_bytes()).expect("coudn't write to output file");

        // reference position left mate
        tag_str = "alnpos_left";
        tag_typeid = 3u8;
        libradicl::rad_types::write_str_bin(tag_str, &libradicl::rad_types::RadIntId::U16, &mut data);
        data.write_all(&tag_typeid.to_le_bytes()).expect("coudn't write to output file");

        // reference position right mate
        tag_str = "alnpos_right";
        tag_typeid = 3u8;
        libradicl::rad_types::write_str_bin(tag_str, &libradicl::rad_types::RadIntId::U16, &mut data);
        data.write_all(&tag_typeid.to_le_bytes()).expect("coudn't write to output file");

        // fragment length estimated based on the alignment
        tag_str = "fraglen";
        tag_typeid = 3u8;
        libradicl::rad_types::write_str_bin(tag_str, &libradicl::rad_types::RadIntId::U16, &mut data);
        data.write_all(&tag_typeid.to_le_bytes()).expect("coudn't write to output file");

        ///////////////////////////////////////// file-level tag values
        // nothing for bulk!
    }

    // dump current buffer content
    owriter.write_all(data.get_ref()).unwrap();

    // let mut bam_reader = bam::Reader::from_path(&input_bam_filename).unwrap();
    // let header_view = bam_reader.header().to_owned();

    // if *threads_count > 1 {
    //     bam_reader.set_threads(threads_count - 1).unwrap();
    // } else {
    //     bam_reader.set_threads(1).unwrap();
    // }

    let reader_threads: Option<usize>;
    if *threads_count > 1 {
        log::info!("thread count: {}", threads_count);
        reader_threads = Some(threads_count - 1);
    } else {
        reader_threads = None;
        log::info!("thread count: {}", threads_count);
    }

    // setup the input BAM Reader
    let mut bqr = BAMQueryRecordReader::new(input_bam_filename, reader_threads);
    let input_header = bqr.get_header().to_owned();

    let mut total_num_chunks = 0u64;
    let mut chunk_reads = 0u32;

    // allocate buffer
    let buf_limit = 10000u32;
    data = Cursor::new(Vec::<u8>::with_capacity((buf_limit * 100) as usize));
    // placeholder for number of bytes and number of records
    data.write_all(&chunk_reads.to_le_bytes()).unwrap();
    data.write_all(&chunk_reads.to_le_bytes()).unwrap();

    // let mut last_qname = String::from("");
    let mut all_query_records: Vec<record::Record> = Vec::new();
    let required_tags: Vec<&str> = vec![];
    // let mut rec_num_align: u32;
    // let mut read_num_align: u32 = 0;
    // let mut first_record: record::Record = record::Record::new();
    // let mut second_record: record::Record; // = record::Record::new();
    // let mut n = 0;

    while let Ok(Some(ret_vec)) = bqr.get_next_query_records() {
        all_query_records.clear();
        for r in ret_vec.iter() {
            let mut txp_records = convert::convert_query_bam_records(r, &input_header, transcripts, txp_lengths, trees, max_softlen, &required_tags);
            all_query_records.append(&mut txp_records);
        }
        if !all_query_records.is_empty() {
            let write_success = dump_collected_alignments_bulk_pe(&all_query_records, &mut data);
            if write_success {
                chunk_reads += 1;
            }
        }
        if chunk_reads >= buf_limit {
            // dump current chunk and start a new one
            data.set_position(0);
            // number of bytes
            data.write_all(&(data.get_ref().len() as u32).to_le_bytes()).unwrap();
            // number reads
            data.write_all(&chunk_reads.to_le_bytes()).unwrap();
            owriter.write_all(data.get_ref()).unwrap();
            total_num_chunks += 1;
            chunk_reads = 0;
            data = Cursor::new(Vec::<u8>::with_capacity((buf_limit * 100) as usize));
            // placeholder for number of bytes and number of records
            data.write_all(&chunk_reads.to_le_bytes()).unwrap();
            data.write_all(&chunk_reads.to_le_bytes()).unwrap();
        }
    }
    // dump the last chunk
    if chunk_reads > 0 {
        data.set_position(0);
        // number of bytes
        data.write_all(&(data.get_ref().len() as u32).to_le_bytes()).unwrap();
        // number reads
        data.write_all(&chunk_reads.to_le_bytes()).unwrap();
        owriter.write_all(data.get_ref()).unwrap();
        total_num_chunks += 1;
    }
    // update the number of chunks in the header
    owriter.flush().expect("File buffer could not be flushed");
    owriter.seek(SeekFrom::Start(end_header_pos)).expect("couldn't seek in output file");
    owriter
        .write_all(&total_num_chunks.to_le_bytes())
        .expect("couldn't write to output file.");
}

fn dump_collected_alignments_singlecell(
    all_read_records: &[record::Record],
    bc_typeid: &u8,
    umi_typeid: &u8,
    corrected_tags: bool,
    owriter: &mut Cursor<Vec<u8>>,
) -> bool {
    // add stored data to the current chunk
    let mut wrote_some: bool = false;
    let mut written_records: u32 = 0;
    let mut data = Cursor::new(vec![]);

    // check if barcode and UMI can be converted to numbers
    let bc_string: String;
    let umi_string: String;
    if corrected_tags {
        if let Ok(Aux::String(bc_str)) = all_read_records.first().unwrap().aux(b"CB") {
            bc_string = bc_str.to_string();
        } else {
            panic!("Input record missing CB tag!");
        }

        if let Ok(Aux::String(umi_str)) = all_read_records.first().unwrap().aux(b"UB") {
            umi_string = umi_str.to_string();
        } else {
            panic!("Input record missing UB tag!");
        }
    } else {
        if let Ok(Aux::String(bc_str)) = all_read_records.first().unwrap().aux(b"CR") {
            bc_string = bc_str.to_string();
        } else {
            panic!("Input record missing CR tag!");
        }

        if let Ok(Aux::String(umi_str)) = all_read_records.first().unwrap().aux(b"UR") {
            umi_string = umi_str.to_string();
        } else {
            panic!("Input record missing UR tag!");
        }
    }

    log::debug!(
        "qname:{} bc:{} umi:{}",
        String::from_utf8(all_read_records.first().unwrap().qname().to_vec()).unwrap(),
        bc_string,
        umi_string
    );
    if (*bc_typeid != 8 && bc_string.contains('N')) || (*umi_typeid != 8 && umi_string.contains('N')) {
        log::debug!("barcode or UMI has N");
        return false;
    }

    for txp_rec in all_read_records.iter() {
        if !txp_rec.is_unmapped() {
            // alignment reference ID
            // data.write_all(&(txp_rec.tid() as u32).to_le_bytes()).unwrap();
            // alignment orientation
            // data.write_all(&(txp_rec.is_reverse() as u8).to_le_bytes()).unwrap();
            // alignment position
            // data.write_all(&(txp_rec.pos() as u32).to_le_bytes()).unwrap();

            // array of alignment-specific tags
            let mut tid_comressed = txp_rec.tid() as u32;
            if !txp_rec.is_reverse() {
                tid_comressed |= 0x80000000_u32;
            }
            data.write_all(&tid_comressed.to_le_bytes()).unwrap();
            written_records += 1;
            wrote_some = true;
        }
    }

    if written_records > 0 {
        // add stored data to the current chunk
        // number of alignments
        owriter.write_all(&written_records.to_le_bytes()).unwrap();
        // read length
        // owriter.write_all(&(all_read_records.first().unwrap().seq_len() as u16).to_le_bytes()).unwrap();
        // read-level tags for single-cell
        // bc
        if *bc_typeid == 8 {
            // write as a string
            libradicl::rad_types::write_str_bin(&bc_string, &libradicl::rad_types::RadIntId::U16, owriter);
        } else {
            // convert to integer
            let bc_int: u64 = cb_string_to_u64(bc_string.as_bytes()).unwrap();
            match decode_int_type_tag(*bc_typeid).unwrap() {
                libradicl::rad_types::RadIntId::U32 => {
                    owriter.write_all(&(bc_int as u32).to_le_bytes()).unwrap();
                }
                libradicl::rad_types::RadIntId::U64 => {
                    owriter.write_all(&(bc_int as u64).to_le_bytes()).unwrap();
                }
                _ => {} // bc_typeid can only be 3, 4, or 8
            }
        }

        // umi
        if *umi_typeid == 8 {
            // write as a string
            libradicl::rad_types::write_str_bin(&umi_string, &libradicl::rad_types::RadIntId::U16, owriter);
        } else {
            // convert to integer
            let umi_int: u64 = cb_string_to_u64(umi_string.as_bytes()).unwrap();
            match decode_int_type_tag(*umi_typeid).unwrap() {
                libradicl::rad_types::RadIntId::U32 => {
                    owriter.write_all(&(umi_int as u32).to_le_bytes()).unwrap();
                }
                libradicl::rad_types::RadIntId::U64 => {
                    owriter.write_all(&(umi_int as u64).to_le_bytes()).unwrap();
                }
                _ => {} // bc_typeid can only be 3, 4, or 8
            }
        }

        // dump records
        owriter.write_all(data.get_ref()).unwrap();
    }

    wrote_some
}

pub fn bam2rad_singlecell(
    input_bam_filename: &String,
    output_dirname: &String,
    rad_mapped_filename: &String,
    rad_unmapped_filename: &String,
    transcripts: &Vec<String>,
    txp_lengths: &Vec<i32>,
    trees: &FnvHashMap<String, COITree<ExonNode, u32>>,
    threads_count: &usize,
    max_softlen: &usize,
    corrected_tags: bool,
) {
    let out_dir_path = Path::new(output_dirname);
    fs::create_dir_all(out_dir_path).unwrap();
    let out_rad_path = out_dir_path.join(rad_mapped_filename);
    let ofile = File::create(out_rad_path.to_str().unwrap()).unwrap();
    let out_unmapped_path = out_dir_path.join(rad_unmapped_filename);
    let _ofile_unmapped = File::create(out_unmapped_path.to_str().unwrap()).unwrap();
    // file writer and intermediate buffer
    let mut owriter = BufWriter::with_capacity(1048576, ofile);
    let mut data = Cursor::new(vec![]);
    //
    let is_paired: u8;
    let end_header_pos: u64;

    let bc_typeid: u8;
    let umi_typeid: u8;

    {
        let first_record = bam_peek(input_bam_filename);
        ////////////////////////////////////////////////// header section
        // NOTE: we assume that all the single-cell protocols are single-end
        // is paired-end
        is_paired = 0;

        data.write_all(&is_paired.to_le_bytes()).expect("couldn't write to output file");
        // number of references
        let ref_count = transcripts.len() as u64;
        log::info!("Number of reference sequences: {}", ref_count);
        data.write_all(&ref_count.to_le_bytes()).expect("couldn't write to output file");
        // list of reference names
        for tname in transcripts.iter() {
            // let name_len = tname.len() as u16;
            // data.write_all(&name_len.to_le_bytes())
            //     .expect("coudn't write to output file");
            // data.write_all(tname.as_bytes())
            //     .expect("coudn't write to output file");
            libradicl::rad_types::write_str_bin(tname, &libradicl::rad_types::RadIntId::U16, &mut data);
        }
        // keep a pointer to header pos
        // end_header_pos = data.seek(SeekFrom::Current(0)).unwrap() - std::mem::size_of::<u64>() as u64;
        end_header_pos = data.seek(SeekFrom::Current(0)).unwrap();
        log::debug!("end header pos: {:?}", end_header_pos);
        // number of chunks
        let initial_num_chunks = 0u64;
        data.write_all(&initial_num_chunks.to_le_bytes()).expect("coudn't write to output file");

        ///////////////////////////////////////// tag definition section
        // FILE-LEVEL tags

        // for single-cell, keep length of cell-barcode and UMI
        let mut num_tags = 2u16;
        data.write_all(&num_tags.to_le_bytes()).expect("coudn't write to output file");

        let mut typeid = 2u8; // u16
        let mut cb_tag_str = "cblen";
        let mut umi_tag_str = "ulen";

        // str - type
        libradicl::rad_types::write_str_bin(cb_tag_str, &libradicl::rad_types::RadIntId::U16, &mut data);
        data.write_all(&typeid.to_le_bytes()).expect("coudn't write to output file");

        // str - type
        libradicl::rad_types::write_str_bin(umi_tag_str, &libradicl::rad_types::RadIntId::U16, &mut data);
        data.write_all(&typeid.to_le_bytes()).expect("coudn't write to output file");

        // READ-LEVEL tags
        let bc_string_in: String;
        let umi_string_in: String;
        if corrected_tags {
            if let Ok(Aux::String(bcs)) = first_record.aux(b"CB") {
                bc_string_in = bcs.to_string();
            } else {
                panic!("Input record missing CB tag!");
            }
            if let Ok(Aux::String(umis)) = first_record.aux(b"UB") {
                umi_string_in = umis.to_string();
            } else {
                panic!("Input record missing UB tag!");
            }
        } else {
            if let Ok(Aux::String(bcs)) = first_record.aux(b"CR") {
                bc_string_in = bcs.to_string();
            } else {
                panic!("Input record missing CR tag!");
            }
            if let Ok(Aux::String(umis)) = first_record.aux(b"UR") {
                umi_string_in = umis.to_string();
            } else {
                panic!("Input record missing UR tag!");
            }
        }

        let bclen = bc_string_in.len() as u16;
        let umilen = umi_string_in.len() as u16;

        num_tags = 2u16;
        data.write_all(&num_tags.to_le_bytes()).expect("coudn't write to output file");
        cb_tag_str = "b";
        umi_tag_str = "u";

        // type is conditional on barcode and umi length
        // this follows SalmonAlevin.cpp implementation
        bc_typeid = match bclen {
            1..=16 => libradicl::rad_types::encode_type_tag(libradicl::rad_types::RadType::U32).unwrap(),
            17..=32 => libradicl::rad_types::encode_type_tag(libradicl::rad_types::RadType::U64).unwrap(),
            _ => 8,
        };

        umi_typeid = match umilen {
            1..=16 => libradicl::rad_types::encode_type_tag(libradicl::rad_types::RadType::U32).unwrap(),
            17..=32 => libradicl::rad_types::encode_type_tag(libradicl::rad_types::RadType::U64).unwrap(),
            _ => 8,
        };

        log::debug!("CB LEN : {}, CB TYPE : {}", bclen, bc_typeid);
        log::debug!("UMI LEN : {}, UMI TYPE : {}", umilen, umi_typeid);

        libradicl::rad_types::write_str_bin(cb_tag_str, &libradicl::rad_types::RadIntId::U16, &mut data);
        data.write_all(&bc_typeid.to_le_bytes()).expect("coudn't write to output file");

        libradicl::rad_types::write_str_bin(umi_tag_str, &libradicl::rad_types::RadIntId::U16, &mut data);
        data.write_all(&umi_typeid.to_le_bytes()).expect("coudn't write to output file");

        // ALIGNMENT-LEVEL tags
        num_tags = 1u16;
        data.write_all(&num_tags.to_le_bytes()).expect("couldn't write to output file");

        // reference id
        let refid_str = "compressed_ori_refid";
        typeid = 3u8;
        libradicl::rad_types::write_str_bin(refid_str, &libradicl::rad_types::RadIntId::U16, &mut data);
        data.write_all(&typeid.to_le_bytes()).expect("coudn't write to output file");

        ///////////////////////////////////////// file-level tag values
        data.write_all(&bclen.to_le_bytes()).expect("coudn't write to output file");
        data.write_all(&umilen.to_le_bytes()).expect("coudn't write to output file");
    }

    // dump current buffer content
    owriter.write_all(data.get_ref()).unwrap();

    // let mut bam_reader = bam::Reader::from_path(&input_bam_filename).unwrap();
    // let header_view = bam_reader.header().to_owned();

    // if *threads_count > 1 {
    //     bam_reader.set_threads(threads_count - 1).unwrap();
    // } else {
    //     bam_reader.set_threads(1).unwrap();
    // }

    let reader_threads: Option<usize>;
    if *threads_count > 1 {
        log::info!("thread count: {}", threads_count);
        reader_threads = Some(threads_count - 1);
    } else {
        reader_threads = None;
        log::info!("thread count: {}", threads_count);
    }

    // setup the input BAM Reader
    let mut bqr = BAMQueryRecordReader::new(input_bam_filename, reader_threads);
    let input_header = bqr.get_header().to_owned();

    let mut total_num_chunks = 0u64;
    let mut chunk_reads = 0u32;

    // allocate buffer
    let buf_limit = 10000u32;
    data = Cursor::new(Vec::<u8>::with_capacity((buf_limit * 100) as usize));
    // placeholder for number of bytes and number of records
    data.write_all(&chunk_reads.to_le_bytes()).unwrap();
    data.write_all(&chunk_reads.to_le_bytes()).unwrap();

    // let mut last_qname = String::from("");
    let mut all_query_records: Vec<record::Record> = Vec::new();
    let required_tags: Vec<&str> = vec![];
    // let mut n = 0;

    while let Ok(Some(ret_vec)) = bqr.get_next_query_records() {
        all_query_records.clear();
        for r in ret_vec.iter() {
            let mut txp_records = convert::convert_query_bam_records(r, &input_header, transcripts, txp_lengths, trees, max_softlen, &required_tags);
            all_query_records.append(&mut txp_records);
        }
        if !all_query_records.is_empty() {
            let write_success = dump_collected_alignments_singlecell(&all_query_records, &bc_typeid, &umi_typeid, corrected_tags, &mut data);
            if write_success {
                chunk_reads += 1;
            }
        }
        if chunk_reads >= buf_limit {
            // dump current chunk and start a new one
            data.set_position(0);
            // number of bytes
            data.write_all(&(data.get_ref().len() as u32).to_le_bytes()).unwrap();
            // number reads
            data.write_all(&chunk_reads.to_le_bytes()).unwrap();
            owriter.write_all(data.get_ref()).unwrap();
            total_num_chunks += 1;
            chunk_reads = 0;
            data = Cursor::new(Vec::<u8>::with_capacity((buf_limit * 100) as usize));
            // placeholder for number of bytes and number of records
            data.write_all(&chunk_reads.to_le_bytes()).unwrap();
            data.write_all(&chunk_reads.to_le_bytes()).unwrap();
        }
    }
    // dump the last chunk
    if chunk_reads > 0 {
        // dump current chunk and start a new one
        data.set_position(0);
        // number of bytes
        data.write_all(&(data.get_ref().len() as u32).to_le_bytes()).unwrap();
        // number reads
        data.write_all(&chunk_reads.to_le_bytes()).unwrap();
        owriter.write_all(data.get_ref()).unwrap();
        total_num_chunks += 1;
    }
    // update the number of chunks in the header
    owriter.flush().expect("File buffer could not be flushed");
    owriter.seek(SeekFrom::Start(end_header_pos)).expect("couldn't seek in output file");
    owriter
        .write_all(&total_num_chunks.to_le_bytes())
        .expect("couldn't write to output file.");
}
