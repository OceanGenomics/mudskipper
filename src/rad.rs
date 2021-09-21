// use std::io::{stdout, BufReader, BufWriter, Cursor, Seek, SeekFrom, Write};

use crate::annotation;
use crate::convert;

extern crate fnv;

// use std::fs;
use std::fs::File;
// use std::path::Path;
use log::{info, error, debug};

use annotation::ExonNode;

use std::io::{BufWriter, Cursor, Write, Seek, SeekFrom};
use rust_htslib::{bam, bam::Read, bam::record};
use fnv::FnvHashMap;
use coitrees::{COITree};

// struct RadRecord {
//     pub fw: bool,
//     pub qname: String,
//     pub tid_list: Vec<i32>
// }

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

pub fn bam2rad(input_bam_filename: &String, 
                 output_rad_filename: &String,
                 transcripts: &Vec<String>,
                 txp_lengths: &Vec<i32>,
                 trees: &FnvHashMap::<String, COITree<ExonNode, u32>>,
                 threads_count: &usize,
                 max_softlen: &usize) {
    
    let ofile = File::create(&output_rad_filename).unwrap();
    // {
    //     let bam_bytes = fs::metadata(&input_bam_filename).unwrap().len();
    //     info!("Bam file size in bytes {:?}", bam_bytes);
    // }
    // file writer and intermediate buffer
    let mut owriter = BufWriter::with_capacity(1048576, ofile);
    let mut data = Cursor::new(vec![]);
    // 
    let is_paired: u8;
    let end_header_pos: u64;

    {
        let first_rec = bam_peek(input_bam_filename);
        ////////////////////////////////////////////////// header section
        // is paired-end
        if first_rec.is_paired() {
            is_paired = 1;
        } else {
            is_paired = 0;
        }

        data.write_all(&is_paired.to_le_bytes())
            .expect("couldn't write to output file");
        // number of references
        let ref_count = transcripts.len() as u64;
        println!("ref_count: {}", ref_count);
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
        info!("end header pos: {:?}", end_header_pos);
        // number of chunks
        let initial_num_chunks = 0u64;
        data.write_all(&initial_num_chunks.to_le_bytes())
            .expect("coudn't write to output file");

        ///////////////////////////////////////// tag definition section
        // FILE-LEVEL tags
        // TODO: assumes bulk sample only for now
        let mut num_tags = 0u16; // no file-level tags for bulk
        data.write_all(&num_tags.to_le_bytes())
            .expect("coudn't write to output file");

        // READ-LEVEL tags
        // TODO: assumes bulk sample only for now
        num_tags = 0u16; // no read-level tags for bulk
        data.write_all(&num_tags.to_le_bytes())
            .expect("coudn't write to output file");

        // ALIGNMENT-LEVEL tags
        // TODO: assumes bulk sample only for now
        num_tags = 0u16;
        data.write_all(&num_tags.to_le_bytes())
            .expect("couldn't write to output file");

        ///////////////////////////////////////// file-level tag values
        // nothing for bulk!
    }

    // TODO: remove dummy bytes
    data.write_all(&(0xffffffffffffffff as u64).to_le_bytes()).expect("coudn't write to output file");

    // dump current buffer content
    owriter.write_all(data.get_ref()).unwrap();

    let mut bam_reader = bam::Reader::from_path(&input_bam_filename).unwrap();
    let header_view = bam_reader.header().to_owned();

    if *threads_count > 1 {
        bam_reader.set_threads(threads_count - 1).unwrap();
    } else {
        bam_reader.set_threads(1).unwrap();
    }

    // NOT NEEDED ANYMORE
    // // get the first record for creating flags
    // let mut rec = bam::Record::new();
    // let first_record_exists = bam.read(&mut rec).is_some();
    // if !first_record_exists {
    //     error!("bam file had no records!");
    //     std::process::exit(1);
    // }

    let mut total_num_chunks = 0u64;
    let mut chunk_reads = 0u32;

    // allocate buffer
    let buf_limit = 10000u32;
    data = Cursor::new(Vec::<u8>::with_capacity((buf_limit * 100) as usize));
    // placeholder for number of bytes and number of records
    data.write_all(&chunk_reads.to_le_bytes()).unwrap();
    data.write_all(&chunk_reads.to_le_bytes()).unwrap();
    
    // TODO: remove dummy bytes
    data.write_all(&(0xffffffffffffffff as u64).to_le_bytes()).expect("coudn't write to output file");

    if is_paired == 0 {
        println!("debug 1: single-end");
        let mut last_qname = String::from("");
        let mut all_read_records: Vec<record::Record> = Vec::new();
        let mut n = 0;
        for rec in bam_reader.records() {
            n = n + 1;
            let record = rec.unwrap();
            let qname = String::from_utf8(record.qname().to_vec()).unwrap();
            debug!("qname: {}", qname);
            println!("qname: {}", qname);
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
                    println!("### dumping for {}, size: {}", last_qname, all_read_records.len());
                    // TODO: remove dummy bytes
                    data.write_all(&(0xaaaaaaaaaaaaaaaa as u64).to_le_bytes()).expect("coudn't write to output file");
                    // add stored data to the current chunk
                    // number of alignments
                    data.write_all(&(all_read_records.len() as u32).to_le_bytes()).unwrap();
                    // read length
                    data.write_all(&(all_read_records.first().unwrap().seq_len() as u16).to_le_bytes()).unwrap();
                    // No read-level tags for bulk

                    for txp_rec in all_read_records.iter() {
                        // alignment reference ID
                        data.write_all(&(txp_rec.tid() as u32).to_le_bytes()).unwrap();
                        // alignment orientation
                        data.write_all(&(txp_rec.is_reverse() as u8).to_le_bytes()).unwrap();
                        // alignment position
                        data.write_all(&(txp_rec.pos() as u32).to_le_bytes()).unwrap();
                        // array of alignment-specific tags
                    }
                    
                    chunk_reads += 1;
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
                    
                    // TODO: remove dummy bytes
                    data.write_all(&(0xffffffffffffffff as u64).to_le_bytes()).expect("coudn't write to output file");
                }
            }
        }
        // add stored data to the last chunk
        if all_read_records.len() > 0 {
            println!("### dumping for {}, size: {}", last_qname, all_read_records.len());
            // TODO: remove dummy bytes
            data.write_all(&(0xaaaaaaaaaaaaaaaa as u64).to_le_bytes()).expect("coudn't write to output file");
            // add stored data to the current chunk
            // number of alignments
            data.write_all(&(all_read_records.len() as u32).to_le_bytes()).unwrap();
            // read length
            data.write_all(&(all_read_records.first().unwrap().seq_len() as u16).to_le_bytes()).unwrap();
            // No read-level tags for bulk

            for txp_rec in all_read_records.iter() {
                // alignment reference ID
                data.write_all(&(txp_rec.tid() as u32).to_le_bytes()).unwrap();
                // alignment orientation
                data.write_all(&(txp_rec.is_reverse() as u8).to_le_bytes()).unwrap();
                // alignment position
                data.write_all(&(txp_rec.pos() as u32).to_le_bytes()).unwrap();
                // array of alignment-specific tags
            }
            
            chunk_reads += 1;
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

    } else {
        println!("debug 1: paired-end");
        let mut last_qname = String::from("");
        let mut all_read_records: Vec<record::Record> = Vec::new();
        let mut rec_num_align: u32;
        let mut read_num_align: u32 = 0;
        let mut first_record: record::Record = record::Record::new();
        let mut second_record: record::Record = record::Record::new();
        let mut n = 0;
        for rec in bam_reader.records() {
            n = n + 1;
            let record = rec.unwrap();
            let qname = String::from_utf8(record.qname().to_vec()).unwrap();
            debug!("qname: {}", qname);
            println!("qname: {}", qname);

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
                println!("rec_num_align: {}", txp_records.len());
                rec_num_align = txp_records.len() as u32;
            } else if second_record.is_unmapped() {
                txp_records = convert::convert_single_end(&first_record,
                                                        &header_view,
                                                        transcripts,
                                                        trees,
                                                        max_softlen);
                println!("rec_num_align: {}", txp_records.len());
                rec_num_align = txp_records.len() as u32;
            } else { // both mates in pair are mapped
                txp_records = convert::convert_paired_end(&first_record,
                                                                    &second_record,
                                                                    &header_view,
                                                                    transcripts,
                                                                    txp_lengths,
                                                                    trees,
                                                                    max_softlen);
                println!("rec_num_align: {}", (txp_records.len() / 2));
                rec_num_align = (txp_records.len() / 2) as u32;
            };
            // num_align += 1;

            if qname == last_qname {
                all_read_records.append(&mut txp_records);
                read_num_align += rec_num_align;
            } else {
                if all_read_records.len() > 0 {
                    println!("### dumping for {}, size: {}, read_num_align: {}", last_qname, all_read_records.len(), read_num_align);
                    // TODO: remove dummy bytes
                    data.write_all(&(0xaaaaaaaaaaaaaaaa as u64).to_le_bytes()).expect("coudn't write to output file");
                    // add stored data to the current chunk
                    // number of alignments
                    data.write_all(&read_num_align.to_le_bytes()).unwrap();
                    // first mate length
                    data.write_all(&(first_record.seq_len() as u16).to_le_bytes()).unwrap();
                    // second mate length
                    data.write_all(&(second_record.seq_len() as u16).to_le_bytes()).unwrap();
                    // No read-level tags for bulk
                    
                    let mut rec_iter = all_read_records.iter();
                    while let Some(rec1) = rec_iter.next() {
                        if !rec1.is_mate_unmapped() { // there should be another mate
                            if let Some(rec2) = rec_iter.next() {
                                // alignment reference ID
                                assert_eq!(rec1.tid(), rec2.tid());
                                data.write_all(&(rec1.tid() as u32).to_le_bytes()).unwrap();
                                let mut aln_type: u8 = 0;
                                let pos_left: u32;
                                let pos_right: u32;
                                if rec1.is_first_in_template() == true { // rec1 is left
                                    if rec1.is_reverse() {
                                        aln_type = aln_type | (2 as u8);
                                    }
                                    if rec2.is_reverse() {
                                        aln_type = aln_type | (1 as u8);
                                    }
                                    pos_left = rec1.pos() as u32;
                                    pos_right = rec2.pos() as u32;
                                } else { // rec2 is left
                                    if rec1.is_reverse() {
                                        aln_type = aln_type | (1 as u8);
                                    }
                                    if rec2.is_reverse() {
                                        aln_type = aln_type | (2 as u8);
                                    }
                                    pos_right = rec1.pos() as u32;
                                    pos_left = rec2.pos() as u32;
                                }
                                // alignment type (0..7)
                                data.write_all(&aln_type.to_le_bytes()).unwrap();
                                // alignment position left mate
                                data.write_all(&pos_left.to_le_bytes()).unwrap();
                                // alignment position right mate
                                data.write_all(&pos_right.to_le_bytes()).unwrap();
                                // array of alignment-specific tags
                            } else {
                                error!("couldn't find respective mate!");
                            }
                        } else { // there is no mate
                            // alignment reference ID
                            data.write_all(&(rec1.tid() as u32).to_le_bytes()).unwrap();
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
                            data.write_all(&aln_type.to_le_bytes()).unwrap();
                            // alignment position of left/right mate
                            data.write_all(&(rec1.pos() as u32).to_le_bytes()).unwrap();
                            // array of alignment-specific tags
                        }
                    }
                    chunk_reads += 1;
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
                    
                    // TODO: remove dummy bytes
                    data.write_all(&(0xffffffffffffffff as u64).to_le_bytes()).expect("coudn't write to output file");
                }
            }
        }
        if all_read_records.len() > 0 {
            println!("### dumping for {}, size: {}, read_num_align: {}", last_qname, all_read_records.len(), read_num_align);
            // TODO: remove dummy bytes
            data.write_all(&(0xaaaaaaaaaaaaaaaa as u64).to_le_bytes()).expect("coudn't write to output file");
            // add stored data to the current chunk
            // number of alignments
            data.write_all(&read_num_align.to_le_bytes()).unwrap();
            // first mate length
            data.write_all(&(first_record.seq_len() as u16).to_le_bytes()).unwrap();
            // second mate length
            data.write_all(&(second_record.seq_len() as u16).to_le_bytes()).unwrap();
            // No read-level tags for bulk
            
            let mut rec_iter = all_read_records.iter();
            while let Some(rec1) = rec_iter.next() {
                if !rec1.is_mate_unmapped() { // there should be another mate
                    if let Some(rec2) = rec_iter.next() {
                        // alignment reference ID
                        assert_eq!(rec1.tid(), rec2.tid());
                        data.write_all(&(rec1.tid() as u32).to_le_bytes()).unwrap();
                        let mut aln_type: u8 = 0;
                        let pos_left: u32;
                        let pos_right: u32;
                        if rec1.is_first_in_template() == true { // rec1 is left
                            if rec1.is_reverse() {
                                aln_type = aln_type | (2 as u8);
                            }
                            if rec2.is_reverse() {
                                aln_type = aln_type | (1 as u8);
                            }
                            pos_left = rec1.pos() as u32;
                            pos_right = rec2.pos() as u32;
                        } else { // rec2 is left
                            if rec1.is_reverse() {
                                aln_type = aln_type | (1 as u8);
                            }
                            if rec2.is_reverse() {
                                aln_type = aln_type | (2 as u8);
                            }
                            pos_right = rec1.pos() as u32;
                            pos_left = rec2.pos() as u32;
                        }
                        // alignment type (0..7)
                        data.write_all(&aln_type.to_le_bytes()).unwrap();
                        // alignment position left mate
                        data.write_all(&pos_left.to_le_bytes()).unwrap();
                        // alignment position right mate
                        data.write_all(&pos_right.to_le_bytes()).unwrap();
                        // array of alignment-specific tags
                    } else {
                        error!("couldn't find respective mate!");
                    }
                } else { // there is no mate
                    // alignment reference ID
                    data.write_all(&(rec1.tid() as u32).to_le_bytes()).unwrap();
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
                    data.write_all(&aln_type.to_le_bytes()).unwrap();
                    // alignment position of left/right mate
                    data.write_all(&(rec1.pos() as u32).to_le_bytes()).unwrap();
                    // array of alignment-specific tags
                }
            }
            chunk_reads += 1;
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
}
