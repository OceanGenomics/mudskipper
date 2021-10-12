use crate::annotation;
use crate::convert;

extern crate fnv;

use std::fs::File;
use std::error::Error;
use libradicl::decode_int_type_tag;
use log::{info, error, debug};

use annotation::ExonNode;

use std::io::{BufWriter, Cursor, Write, Seek, SeekFrom};
use rust_htslib::{bam, bam::Read, bam::record, bam::record::Aux};
use fnv::FnvHashMap;
use coitrees::{COITree};
// use scroll::Pread;

// fn read_into_u64<T: std::io::Read>(reader: &mut T, rt: &libradicl::RadIntId) -> u64 {
//     let mut rbuf = [0u8; 8];
//     let v: u64;
//     match rt {
//         libradicl::RadIntId::U8 => {
//             reader.read_exact(&mut rbuf[0..1]).unwrap();
//             v = rbuf.pread::<u8>(0).unwrap() as u64;
//         }
//         libradicl::RadIntId::U16 => {
//             reader.read_exact(&mut rbuf[0..2]).unwrap();
//             v = rbuf.pread::<u16>(0).unwrap() as u64;
//         }
//         libradicl::RadIntId::U32 => {
//             reader.read_exact(&mut rbuf[0..4]).unwrap();
//             v = rbuf.pread::<u32>(0).unwrap() as u64;
//         }
//         libradicl::RadIntId::U64 => {
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
        error!("bam file had no records!");
        std::process::exit(1);
    }
    rec
}

pub fn bam2rad_bulk_wrapper(input_bam_filename: &String,
                 output_rad_filename: &String,
                 transcripts: &Vec<String>,
                 txp_lengths: &Vec<i32>,
                 trees: &FnvHashMap::<String, COITree<ExonNode, u32>>,
                 threads_count: &usize,
                 max_softlen: &usize) {
    let first_record = bam_peek(input_bam_filename);
    ////////////////////////////////////////////////// header section
    // is paired-end
    if first_record.is_paired() {
        bam2rad_bulk_pe(input_bam_filename, output_rad_filename, transcripts, txp_lengths, trees, threads_count, max_softlen);
    } else {
        bam2rad_bulk_se(input_bam_filename, output_rad_filename, transcripts, trees, threads_count, max_softlen);
    }
}

fn dump_collected_alignments_bulk_se(all_read_records: &Vec<record::Record>, owriter: &mut Cursor<Vec<u8>>) -> bool {
    // add stored data to the current chunk
    // number of alignments
    owriter.write_all(&(all_read_records.len() as u32).to_le_bytes()).unwrap();
    // read length
    // owriter.write_all(&(all_read_records.first().unwrap().seq_len() as u16).to_le_bytes()).unwrap();
    // No read-level tags for bulk

    for txp_rec in all_read_records.iter() {
        // alignment reference ID
        // owriter.write_all(&(txp_rec.tid() as u32).to_le_bytes()).unwrap();
        // alignment orientation
        // owriter.write_all(&(txp_rec.is_reverse() as u8).to_le_bytes()).unwrap();
        // alignment position
        // owriter.write_all(&(txp_rec.pos() as u32).to_le_bytes()).unwrap();
        
        // array of alignment-specific tags
        let mut tid_comressed = txp_rec.tid() as u32;
        if txp_rec.is_reverse() == false {
            tid_comressed |= 0x80000000 as u32;
        }
        owriter.write_all(&tid_comressed.to_le_bytes()).unwrap();
    }
    return true;
}

pub fn bam2rad_bulk_se(input_bam_filename: &String,
                 output_rad_filename: &String,
                 transcripts: &Vec<String>,
                 trees: &FnvHashMap::<String, COITree<ExonNode, u32>>,
                 threads_count: &usize,
                 max_softlen: &usize) {
    
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

        data.write_all(&is_paired.to_le_bytes())
            .expect("couldn't write to output file");
        // number of references
        let ref_count = transcripts.len() as u64;
        info!("Number of reference sequences: {}", ref_count);
        data.write_all(&ref_count.to_le_bytes())
            .expect("couldn't write to output file");
        // list of reference names
        for tname in transcripts.iter() {            
            let name_len = tname.len() as u16;
            data.write_all(&name_len.to_le_bytes())
                .expect("coudn't write to output file");
            data.write_all(tname.as_bytes())
                .expect("coudn't write to output file");
        }
        // keep a pointer to header pos
        // end_header_pos = data.seek(SeekFrom::Current(0)).unwrap() - std::mem::size_of::<u64>() as u64;
        end_header_pos = data.seek(SeekFrom::Current(0)).unwrap();
        debug!("end header pos: {:?}", end_header_pos);
        // number of chunks
        let initial_num_chunks = 0u64;
        data.write_all(&initial_num_chunks.to_le_bytes())
            .expect("coudn't write to output file");

        ///////////////////////////////////////// tag definition section
        // FILE-LEVEL tags
        let mut num_tags = 0u16; // no file-level tags for bulk
        data.write_all(&num_tags.to_le_bytes())
            .expect("coudn't write to output file");

        // READ-LEVEL tags
        num_tags = 0u16; // no read-level tags for bulk
        data.write_all(&num_tags.to_le_bytes())
            .expect("coudn't write to output file");

        // ALIGNMENT-LEVEL tags
        num_tags = 1u16;
        data.write_all(&num_tags.to_le_bytes())
            .expect("couldn't write to output file");

        // reference id
        let refid_str = "compressed_ori_refid";
        let typeid = 3u8;
        libradicl::write_str_bin(&refid_str, &libradicl::RadIntId::U16, &mut data);
        data.write_all(&typeid.to_le_bytes())
            .expect("coudn't write to output file");

        ///////////////////////////////////////// file-level tag values
        // nothing for bulk!
    }

    // dump current buffer content
    owriter.write_all(data.get_ref()).unwrap();

    let mut bam_reader = bam::Reader::from_path(&input_bam_filename).unwrap();
    let header_view = bam_reader.header().to_owned();

    if *threads_count > 1 {
        bam_reader.set_threads(threads_count - 1).unwrap();
    } else {
        bam_reader.set_threads(1).unwrap();
    }

    let mut total_num_chunks = 0u64;
    let mut chunk_reads = 0u32;

    // allocate buffer
    let buf_limit = 10000u32;
    data = Cursor::new(Vec::<u8>::with_capacity((buf_limit * 100) as usize));
    // placeholder for number of bytes and number of records
    data.write_all(&chunk_reads.to_le_bytes()).unwrap();
    data.write_all(&chunk_reads.to_le_bytes()).unwrap();

    let mut last_qname = String::from("");
    let mut all_read_records: Vec<record::Record> = Vec::new();
    let mut n = 0;
    for rec in bam_reader.records() {
        n = n + 1;
        let record = rec.unwrap();
        let qname = String::from_utf8(record.qname().to_vec()).unwrap();
        debug!("qname: {}", qname);
        if record.is_unmapped() {
            continue;
        }
        // 
        let mut txp_records = convert::convert_single_end(&record, 
                                                          &header_view,
                                                          transcripts,
                                                          trees,
                                                          max_softlen);
        if qname == last_qname {
            all_read_records.append(&mut txp_records);
        } else {
            if all_read_records.len() > 0 {
                let write_success = dump_collected_alignments_bulk_se(&all_read_records, &mut data);
                if write_success == true {
                    chunk_reads += 1;
                }
            }

            // prepare for the new read
            last_qname = qname;
            all_read_records.clear();
            all_read_records.append(&mut txp_records);

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
    }
    // add stored data to the last chunk
    if all_read_records.len() > 0 {
        let write_success = dump_collected_alignments_bulk_se(&all_read_records, &mut data);
        if write_success == true {
            chunk_reads += 1;
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
    owriter.seek(SeekFrom::Start(end_header_pos))
        .expect("couldn't seek in output file");
    owriter.write_all(&total_num_chunks.to_le_bytes())
        .expect("couldn't write to output file.");
}

fn dump_collected_alignments_bulk_pe(all_read_records: &Vec<record::Record>,
                                     read_num_align: u32,
                                     owriter: &mut Cursor<Vec<u8>>) -> bool {
    // add stored data to the current chunk
    // number of alignments
    owriter.write_all(&read_num_align.to_le_bytes()).unwrap();
    // first mate length
    // owriter.write_all(&(first_record.seq_len() as u16).to_le_bytes()).unwrap();
    // second mate length
    // owriter.write_all(&(second_record.seq_len() as u16).to_le_bytes()).unwrap();
    // No read-level tags for bulk
    
    let mut rec_iter = all_read_records.iter();
    while let Some(rec1) = rec_iter.next() {
        if !rec1.is_mate_unmapped() { // there should be another mate
            if let Some(rec2) = rec_iter.next() {
                // alignment reference ID
                assert_eq!(rec1.tid(), rec2.tid());
                owriter.write_all(&(rec1.tid() as u32).to_le_bytes()).unwrap();
                let mut aln_type: u8 = 0;
                // let pos_left: u32;
                // let pos_right: u32;
                if rec1.is_first_in_template() == true { // rec1 is left
                    if rec1.is_reverse() {
                        aln_type = aln_type | (2 as u8);
                    }
                    if rec2.is_reverse() {
                        aln_type = aln_type | (1 as u8);
                    }
                    // pos_left = rec1.pos() as u32;
                    // pos_right = rec2.pos() as u32;
                } else { // rec2 is left
                    if rec1.is_reverse() {
                        aln_type = aln_type | (1 as u8);
                    }
                    if rec2.is_reverse() {
                        aln_type = aln_type | (2 as u8);
                    }
                    // pos_right = rec1.pos() as u32;
                    // pos_left = rec2.pos() as u32;
                }
                // alignment type (0..7)
                owriter.write_all(&aln_type.to_le_bytes()).unwrap();
                // alignment position left mate
                // owriter.write_all(&pos_left.to_le_bytes()).unwrap();
                // alignment position right mate
                // owriter.write_all(&pos_right.to_le_bytes()).unwrap();
                // array of alignment-specific tags
            } else {
                error!("couldn't find respective mate!");
                return false;
            }
        } else { // there is no mate
            // alignment reference ID
            owriter.write_all(&(rec1.tid() as u32).to_le_bytes()).unwrap();
            let aln_type: u8;
            if rec1.is_first_in_template() == true { // right is unmapped
                aln_type = if rec1.is_reverse() {
                    5
                } else {
                    4
                };
            } else { // left is unmapped
                aln_type = if rec1.is_reverse() {
                    7
                } else {
                    6
                };
            }
            // alignment type (0..7)
            owriter.write_all(&aln_type.to_le_bytes()).unwrap();
            // alignment position of left/right mate
            // owriter.write_all(&(rec1.pos() as u32).to_le_bytes()).unwrap();
            // array of alignment-specific tags
        }
    }
    return true;
}

pub fn bam2rad_bulk_pe(input_bam_filename: &String,
                 output_rad_filename: &String,
                 transcripts: &Vec<String>,
                 txp_lengths: &Vec<i32>,
                 trees: &FnvHashMap::<String, COITree<ExonNode, u32>>,
                 threads_count: &usize,
                 max_softlen: &usize) {
    
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

        data.write_all(&is_paired.to_le_bytes())
            .expect("couldn't write to output file");
        // number of references
        let ref_count = transcripts.len() as u64;
        info!("Number of reference sequences: {}", ref_count);
        data.write_all(&ref_count.to_le_bytes())
            .expect("couldn't write to output file");
        // list of reference names
        for tname in transcripts.iter() {            
            let name_len = tname.len() as u16;
            data.write_all(&name_len.to_le_bytes())
                .expect("coudn't write to output file");
            data.write_all(tname.as_bytes())
                .expect("coudn't write to output file");
        }
        // keep a pointer to header pos
        // end_header_pos = data.seek(SeekFrom::Current(0)).unwrap() - std::mem::size_of::<u64>() as u64;
        end_header_pos = data.seek(SeekFrom::Current(0)).unwrap();
        debug!("end header pos: {:?}", end_header_pos);
        // number of chunks
        let initial_num_chunks = 0u64;
        data.write_all(&initial_num_chunks.to_le_bytes())
            .expect("coudn't write to output file");

        ///////////////////////////////////////// tag definition section
        // FILE-LEVEL tags
        let mut num_tags = 0u16; // no file-level tags for bulk
        data.write_all(&num_tags.to_le_bytes())
            .expect("coudn't write to output file");

        // READ-LEVEL tags
        num_tags = 0u16; // no read-level tags for bulk
        data.write_all(&num_tags.to_le_bytes())
            .expect("coudn't write to output file");

        // ALIGNMENT-LEVEL tags
        num_tags = 2u16;
        data.write_all(&num_tags.to_le_bytes())
            .expect("couldn't write to output file");

        // reference id
        let refid_str = "refid";
        let mut typeid = 3u8;
        libradicl::write_str_bin(&refid_str, &libradicl::RadIntId::U16, &mut data);
        data.write_all(&typeid.to_le_bytes())
            .expect("coudn't write to output file");

        // alignment type
        let alntype_str = "alntype";
        typeid = 1u8;
        libradicl::write_str_bin(&alntype_str, &libradicl::RadIntId::U16, &mut data);
        data.write_all(&typeid.to_le_bytes())
            .expect("coudn't write to output file");

        ///////////////////////////////////////// file-level tag values
        // nothing for bulk!
    }

    // dump current buffer content
    owriter.write_all(data.get_ref()).unwrap();

    let mut bam_reader = bam::Reader::from_path(&input_bam_filename).unwrap();
    let header_view = bam_reader.header().to_owned();

    if *threads_count > 1 {
        bam_reader.set_threads(threads_count - 1).unwrap();
    } else {
        bam_reader.set_threads(1).unwrap();
    }

    let mut total_num_chunks = 0u64;
    let mut chunk_reads = 0u32;

    // allocate buffer
    let buf_limit = 10000u32;
    data = Cursor::new(Vec::<u8>::with_capacity((buf_limit * 100) as usize));
    // placeholder for number of bytes and number of records
    data.write_all(&chunk_reads.to_le_bytes()).unwrap();
    data.write_all(&chunk_reads.to_le_bytes()).unwrap();

    let mut last_qname = String::from("");
    let mut all_read_records: Vec<record::Record> = Vec::new();
    let mut rec_num_align: u32;
    let mut read_num_align: u32 = 0;
    let mut first_record: record::Record = record::Record::new();
    let mut second_record: record::Record; // = record::Record::new();
    let mut n = 0;
    for rec in bam_reader.records() {
        n = n + 1;
        let record = rec.unwrap();
        let qname = String::from_utf8(record.qname().to_vec()).unwrap();
        debug!("qname: {}", qname);

        if n % 2 == 1 {
            // this is the first read in pair... save and wait for the second read in the pair
            first_record = record.clone();
            continue;
        }

        // this is the second read in pair... 
        second_record = record.clone();
        assert_eq!(first_record.qname(), second_record.qname());

        let mut txp_records: Vec<record::Record>;
        if first_record.is_unmapped() {
            txp_records = convert::convert_single_end(&second_record,
                                                    &header_view,
                                                    transcripts,
                                                    trees,
                                                    max_softlen);
            // debug!("rec_num_align: {}", txp_records.len());
            rec_num_align = txp_records.len() as u32;
        } else if second_record.is_unmapped() {
            txp_records = convert::convert_single_end(&first_record,
                                                    &header_view,
                                                    transcripts,
                                                    trees,
                                                    max_softlen);
            // debug!("rec_num_align: {}", txp_records.len());
            rec_num_align = txp_records.len() as u32;
        } else { // both mates in pair are mapped
            txp_records = convert::convert_paired_end(&first_record,
                                                                &second_record,
                                                                &header_view,
                                                                transcripts,
                                                                txp_lengths,
                                                                trees,
                                                                max_softlen);
            // debug!("rec_num_align: {}", (txp_records.len() / 2));
            rec_num_align = (txp_records.len() / 2) as u32;
        };

        if qname == last_qname {
            all_read_records.append(&mut txp_records);
            read_num_align += rec_num_align;
        } else {
            if all_read_records.len() > 0 {
                let write_success = dump_collected_alignments_bulk_pe(&all_read_records, read_num_align, &mut data);
                if write_success == true {
                    chunk_reads += 1;
                }
            }

            // prepare for the new read
            last_qname = qname;
            all_read_records.clear();
            all_read_records.append(&mut txp_records);
            read_num_align = rec_num_align;

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
    }
    if all_read_records.len() > 0 {
        let write_success = dump_collected_alignments_bulk_pe(&all_read_records, read_num_align, &mut data);
        if write_success == true {
            chunk_reads += 1;
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
    owriter.seek(SeekFrom::Start(end_header_pos))
        .expect("couldn't seek in output file");
    owriter.write_all(&total_num_chunks.to_le_bytes())
        .expect("couldn't write to output file.");
}

fn dump_collected_alignments_singlecell(all_read_records: &Vec<record::Record>,
                                    bc_typeid: &u8,
                                    umi_typeid: &u8,
                                    owriter: &mut Cursor<Vec<u8>>) -> bool {
    // check if barcode and UMI can be converted to numbers
    let bc_string;
    if let Ok(Aux::String(bc_str)) = all_read_records.first().unwrap().aux(b"CR") {
        bc_string = bc_str.to_string();
    } else {
        panic!("Input record missing CR tag!");
    }

    let umi_string;
    if let Ok(Aux::String(umi_str)) = all_read_records.first().unwrap().aux(b"UR") {
        umi_string = umi_str.to_string();
    } else {
        panic!("Input record missing UR tag!");
    }

    debug!("qname:{} bc:{} umi:{}", String::from_utf8(all_read_records.first().unwrap().qname().to_vec()).unwrap(), bc_string, umi_string);
    if (*bc_typeid != 8 && bc_string.contains('N') == true) || (*umi_typeid != 8 && umi_string.contains('N') == true) {
        debug!("barcode or UMI has N");
        return false;
    }

    // add stored data to the current chunk
    // number of alignments
    owriter.write_all(&(all_read_records.len() as u32).to_le_bytes()).unwrap();
    // read length
    // owriter.write_all(&(all_read_records.first().unwrap().seq_len() as u16).to_le_bytes()).unwrap();
    // read-level tags for single-cell
    // bc
    if *bc_typeid == 8 { // write as a string
        libradicl::write_str_bin(&bc_string, &libradicl::RadIntId::U16, owriter);
    } else { // convert to integer
        let bc_int: u64 = cb_string_to_u64(bc_string.as_bytes()).unwrap();
        match decode_int_type_tag(*bc_typeid).unwrap() {
            libradicl::RadIntId::U32 => {
                owriter.write_all(&(bc_int as u32).to_le_bytes()).unwrap();
            }
            libradicl::RadIntId::U64 => {
                owriter.write_all(&(bc_int as u64).to_le_bytes()).unwrap();
            }
            _ => {} // bc_typeid can only be 3, 4, or 8
        }
    }

    // umi
    if *umi_typeid == 8 { // write as a string
        libradicl::write_str_bin(&umi_string, &libradicl::RadIntId::U16, owriter);
    } else { // convert to integer
        let umi_int: u64 = cb_string_to_u64(umi_string.as_bytes()).unwrap();
        match decode_int_type_tag(*umi_typeid).unwrap() {
            libradicl::RadIntId::U32 => {
                owriter.write_all(&(umi_int as u32).to_le_bytes()).unwrap();
            }
            libradicl::RadIntId::U64 => {
                owriter.write_all(&(umi_int as u64).to_le_bytes()).unwrap();
            }
            _ => {} // bc_typeid can only be 3, 4, or 8
        }
    }

    for txp_rec in all_read_records.iter() {
        // alignment reference ID
        // owriter.write_all(&(txp_rec.tid() as u32).to_le_bytes()).unwrap();
        // alignment orientation
        // owriter.write_all(&(txp_rec.is_reverse() as u8).to_le_bytes()).unwrap();
        // alignment position
        // owriter.write_all(&(txp_rec.pos() as u32).to_le_bytes()).unwrap();
        
        // array of alignment-specific tags
        let mut tid_comressed = txp_rec.tid() as u32;
        if txp_rec.is_reverse() == false {
            tid_comressed |= 0x80000000 as u32;
        }
        owriter.write_all(&tid_comressed.to_le_bytes()).unwrap();
    }
    return true;
}

pub fn bam2rad_singlecell(input_bam_filename: &String, 
                 output_rad_filename: &String,
                 transcripts: &Vec<String>,
                 trees: &FnvHashMap::<String, COITree<ExonNode, u32>>,
                 threads_count: &usize,
                 max_softlen: &usize) {
    
    let ofile = File::create(&output_rad_filename).unwrap();
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

        data.write_all(&is_paired.to_le_bytes())
            .expect("couldn't write to output file");
        // number of references
        let ref_count = transcripts.len() as u64;
        info!("Number of reference sequences: {}", ref_count);
        data.write_all(&ref_count.to_le_bytes())
            .expect("couldn't write to output file");
        // list of reference names
        for tname in transcripts.iter() {
            // let name_len = tname.len() as u16;
            // data.write_all(&name_len.to_le_bytes())
            //     .expect("coudn't write to output file");
            // data.write_all(tname.as_bytes())
            //     .expect("coudn't write to output file");
            libradicl::write_str_bin(tname, &libradicl::RadIntId::U16, &mut data);
        }
        // keep a pointer to header pos
        // end_header_pos = data.seek(SeekFrom::Current(0)).unwrap() - std::mem::size_of::<u64>() as u64;
        end_header_pos = data.seek(SeekFrom::Current(0)).unwrap();
        debug!("end header pos: {:?}", end_header_pos);
        // number of chunks
        let initial_num_chunks = 0u64;
        data.write_all(&initial_num_chunks.to_le_bytes())
            .expect("coudn't write to output file");

        ///////////////////////////////////////// tag definition section
        // FILE-LEVEL tags

        // for single-cell, keep length of cell-barcode and UMI
        let mut num_tags = 2u16;
        data.write_all(&num_tags.to_le_bytes())
            .expect("coudn't write to output file");
        
        let mut typeid = 2u8;  // u16
        let mut cb_tag_str = "cblen";
        let mut umi_tag_str = "ulen";

        // str - type
        libradicl::write_str_bin(&cb_tag_str, &libradicl::RadIntId::U16, &mut data);
        data.write_all(&typeid.to_le_bytes())
            .expect("coudn't write to output file");

        // str - type
        libradicl::write_str_bin(&umi_tag_str, &libradicl::RadIntId::U16, &mut data);
        data.write_all(&typeid.to_le_bytes())
            .expect("coudn't write to output file");

        // READ-LEVEL tags
        let bc_string_in;
        if let Ok(Aux::String(bcs)) = first_record.aux(b"CR") {
            bc_string_in = bcs.to_string();
        } else {
            panic!("Input record missing CR tag!");
        }

        let umi_string_in;
        if let Ok(Aux::String(umis)) = first_record.aux(b"UR") {
            umi_string_in = umis.to_string();
        } else {
            panic!("Input record missing UR tag!");
        }

        let bclen = bc_string_in.len() as u16;
        let umilen = umi_string_in.len() as u16;

        num_tags = 2u16;
        data.write_all(&num_tags.to_le_bytes())
            .expect("coudn't write to output file");
        cb_tag_str = "b";
        umi_tag_str = "u";

        // type is conditional on barcode and umi length
        // this follows SalmonAlevin.cpp implementation
        bc_typeid = match bclen {
            1..=16 => libradicl::encode_type_tag(libradicl::RadType::U32).unwrap(),
            17..=32 => libradicl::encode_type_tag(libradicl::RadType::U64).unwrap(),
            _ => 8,
        };

        umi_typeid = match umilen {
            1..=16 => libradicl::encode_type_tag(libradicl::RadType::U32).unwrap(),
            17..=32 => libradicl::encode_type_tag(libradicl::RadType::U64).unwrap(),
            _ => 8,
        };

        debug!("CB LEN : {}, CB TYPE : {}", bclen, bc_typeid);
        debug!("UMI LEN : {}, UMI TYPE : {}", umilen, umi_typeid);

        libradicl::write_str_bin(&cb_tag_str, &libradicl::RadIntId::U16, &mut data);
        data.write_all(&bc_typeid.to_le_bytes())
            .expect("coudn't write to output file");

        libradicl::write_str_bin(&umi_tag_str, &libradicl::RadIntId::U16, &mut data);
        data.write_all(&umi_typeid.to_le_bytes())
            .expect("coudn't write to output file");

        // ALIGNMENT-LEVEL tags
        num_tags = 1u16;
        data.write_all(&num_tags.to_le_bytes())
            .expect("couldn't write to output file");

        // reference id
        let refid_str = "compressed_ori_refid";
        typeid = 3u8;
        libradicl::write_str_bin(&refid_str, &libradicl::RadIntId::U16, &mut data);
        data.write_all(&typeid.to_le_bytes())
            .expect("coudn't write to output file");

        ///////////////////////////////////////// file-level tag values
        data.write_all(&bclen.to_le_bytes())
            .expect("coudn't write to output file");
        data.write_all(&umilen.to_le_bytes())
            .expect("coudn't write to output file");
    }

    // dump current buffer content
    owriter.write_all(data.get_ref()).unwrap();

    let mut bam_reader = bam::Reader::from_path(&input_bam_filename).unwrap();
    let header_view = bam_reader.header().to_owned();

    if *threads_count > 1 {
        bam_reader.set_threads(threads_count - 1).unwrap();
    } else {
        bam_reader.set_threads(1).unwrap();
    }

    let mut total_num_chunks = 0u64;
    let mut chunk_reads = 0u32;

    // allocate buffer
    let buf_limit = 10000u32;
    data = Cursor::new(Vec::<u8>::with_capacity((buf_limit * 100) as usize));
    // placeholder for number of bytes and number of records
    data.write_all(&chunk_reads.to_le_bytes()).unwrap();
    data.write_all(&chunk_reads.to_le_bytes()).unwrap();

    let mut last_qname = String::from("");
    let mut all_read_records: Vec<record::Record> = Vec::new();
    let mut n = 0;
    for rec in bam_reader.records() {
        n = n + 1;
        let record = rec.unwrap();
        let qname = String::from_utf8(record.qname().to_vec()).unwrap();
        debug!("qname: {}", qname);
        if record.is_unmapped() {
            continue;
        }
        // 
        let mut txp_records = convert::convert_single_end(&record, 
                                                                    &header_view,
                                                                    transcripts,
                                                                    trees,
                                                                    max_softlen);
        if qname == last_qname {
            all_read_records.append(&mut txp_records);
        } else {
            if all_read_records.len() > 0 {
                let write_success = dump_collected_alignments_singlecell(&all_read_records, &bc_typeid, &umi_typeid, &mut data);
                if write_success == true {
                    chunk_reads += 1;
                }
            }

            // prepare for the new read
            last_qname = qname;
            all_read_records.clear();
            all_read_records.append(&mut txp_records);

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
    }
    // add stored data to the last chunk
    if all_read_records.len() > 0 {
        let write_success = dump_collected_alignments_singlecell(&all_read_records, &bc_typeid, &umi_typeid, &mut data);
        if write_success == true {
            chunk_reads += 1;
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
    owriter.seek(SeekFrom::Start(end_header_pos))
        .expect("couldn't seek in output file");
    owriter.write_all(&total_num_chunks.to_le_bytes())
        .expect("couldn't write to output file.");
}
