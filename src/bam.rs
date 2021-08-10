use coitrees::{COITree};

use rust_htslib::bam::{Format, Header, Read, Reader, Writer, header, record};

use crate::annotations;
use annotations::ExonNode;

extern crate bio_types;
use bio_types::strand::Strand;

use crate::intersection;

extern crate fnv;
use fnv::FnvHashMap;

use log::{debug};

pub fn bam2bam(input_bam_filename: &String, 
                output_bam_filename: &String,
                transcripts: &Vec<String>,
                txp_lengths: &Vec<i32>,
                trees: &FnvHashMap::<String, COITree<ExonNode, u32>>,
                threads_count: &usize) -> i32 {
    let mut input_bam = Reader::from_path(input_bam_filename).unwrap();
    // let bam_records = input_bam.records();

    // let header_ = Header::from_template(input_bam.header());
    // let header_text = String::from_utf8(header_.to_bytes().to_owned()).unwrap().;

    let input_bam_header = Reader::from_path(input_bam_filename).unwrap();
    let header_view = input_bam_header.header();
    
    let mut new_header = Header::new();
    // new_header.push_record(header::HeaderRecord::new(b"HD").push_tag(b"VN", &"1.4"));

    let mut counter = 0;
    for tid in transcripts.iter() {
        let txp_len = txp_lengths[counter];
        new_header.push_record(header::HeaderRecord::new(b"SQ").push_tag(b"SN", &tid).push_tag(b"LN", &txp_len));
        counter += 1;
    }

    let mut output_bam = Writer::from_path(output_bam_filename, &new_header, Format::Bam).unwrap();

    if *threads_count >= 2 {
        let threads_count_half = threads_count / 2;
        println!("thread count: {}", threads_count_half);
        input_bam.set_threads(threads_count_half).expect("Failed to set number of BAM reading threads.");
        output_bam.set_threads(threads_count_half).expect("Failed to set number of BAM writing threads.");
    }

    let mut n = 0;
    let mut missed_read = 0;
    let mut first_in_pair = true;
    let mut first_record : record::Record = record::Record::new();
    let mut new_cigar : record::CigarString = record::CigarString(vec![record::Cigar::Match(100)]);
    let mut first_new_cigar : record::CigarString = record::CigarString(vec![record::Cigar::Match(100)]);
    let mut len1 = 0;
    let mut len2 = 0;
    for rec in input_bam.records() {
        n = n + 1;
        let record = rec.unwrap();
        debug!("qname: {}", String::from_utf8(record.qname().to_vec()).unwrap());
        if !record.is_paired() {
            let ranges = intersection::find_ranges_single(&(record.pos() as i32),
                                                          &record.cigar(),
                                                          &mut new_cigar,
                                                          &mut len1);
            let genome_tname = String::from_utf8(header_view.tid2name(record.tid() as u32).to_vec())
                                                .expect("cannot find the tname!");
            if let Some(tree) = trees.get(&genome_tname) {
                let tids = intersection::find_tid(&tree, &ranges);
                if tids.len() > 0 {
                    for (tid, pos_strand) in tids.iter() {
                        debug!("{} {}", tid, transcripts[*tid as usize]);
                        debug!("{}", tid);

                        let mut record_ = record.clone();
                        let mut pos = 0;
                        if pos_strand.1 == Strand::Forward {
                            pos = record.pos() - (pos_strand.0 as i64);
                        } else if pos_strand.1 == Strand::Reverse { 
                            pos = (pos_strand.0 as i64) - record.pos() - len1 as i64;
                            if record.is_reverse() {
                                record_.unset_reverse();
                            } else {
                                record_.set_reverse();
                            }
                        }

                        record_.set(record.qname(),
                                    Some(&new_cigar),
                                    &record.seq().as_bytes(),
                                    record.qual());
                        record_.set_tid(*tid);
                        record_.set_pos(pos);

                        output_bam.write(&record_).unwrap();   
                    }
                } else {
                    missed_read = missed_read + 1;
                }
            } else {
                // log for unannotated splicing junction
            }
            debug!("{} {}", ranges[0].0, ranges[0].1);
            
        } else {
            if first_in_pair {
                first_in_pair = false;
                first_record = record;
            } else {
                first_in_pair = true;
                /* let ranges = intersection::find_tids_paired(&(first_record.pos() as i32), 
                                                              &first_record.cigar(),
                                                              &mut first_new_cigar,
                                                              &(record.pos() as i32), 
                                                              &record.cigar(), 
                                                              &mut new_cigar); */
                let genome_tname = String::from_utf8(header_view.tid2name(record.tid() as u32)
                                                                .to_vec()).expect("cannot find the tname!");
                if let Some(tree) = trees.get(&genome_tname) {
                    let tids = intersection::find_tids_paired(&tree,
                                                                &(first_record.pos() as i32), 
                                                                &first_record.cigar(),
                                                                &mut first_new_cigar,
                                                                &mut len1,
                                                                &(record.pos() as i32), 
                                                                &record.cigar(), 
                                                                &mut new_cigar,
                                                                &mut len2);
                    debug!("{}: {}", first_record.cigar(), first_record.cigar().len());
                    debug!("{}: {}", record.cigar(), record.cigar().len());
                    if tids.len() > 0 {
                        for (tid, pos_strand) in tids.iter() {
                            let mut first_record_ = first_record.clone();
                            let mut second_record_ = record.clone();

                            let mut first_pos = 0;
                            let mut second_pos = 0;
                            if pos_strand.0.1 == Strand::Forward {
                                first_pos = first_record.pos() - (pos_strand.0.0 as i64);
                                debug!("first_pos:{} - pos:{} = {}",
                                        first_record.pos(), pos_strand.0.0, first_pos);
                                second_pos = record.pos() - (pos_strand.1.0 as i64);    
                                debug!("second_pos:{} - pos:{} = {}",
                                    record.pos(), pos_strand.1.0, second_pos);
                            } else if pos_strand.0.1 == Strand::Reverse{ 
                                first_pos = (pos_strand.0.0 as i64) - first_record.pos() - len1 as i64;
                                debug!("pos:{} - first_pos:{} - len1:{} = {}",
                                        pos_strand.0.0, first_record.pos(), len1, first_pos);
                                second_pos = (pos_strand.1.0 as i64) - record.pos() - len2 as i64;
                                debug!("pos:{} - second_pos:{} - len2:{} = {}",
                                        pos_strand.1.0, record.pos(), len2, second_pos);
                                if first_record.is_reverse() {
                                    first_record_.unset_reverse();
                                } else {
                                    first_record_.set_reverse();
                                }
                                if record.is_reverse() {
                                    second_record_.unset_reverse();
                                } else {
                                    second_record_.set_reverse();
                                }

                            }
                            debug!("{} {}", tid, transcripts[*tid as usize]);
                            debug!("{}", tid);

                            first_record_.set(first_record.qname(), 
                                              Some(&first_new_cigar), 
                                              &first_record.seq().as_bytes(), 
                                              first_record.qual());
                            first_record_.set_tid(*tid);

                            first_record_.set_pos(first_pos);
                            first_record_.set_mpos(second_pos);

                            output_bam.write(&first_record_).unwrap();

                            second_record_.set(record.qname(), Some(&new_cigar), 
                                               &record.seq().as_bytes(), 
                                               record.qual());
                            second_record_.set_tid(*tid);
                            second_record_.set_pos(second_pos);
                            second_record_.set_mpos(first_pos);

                            output_bam.write(&second_record_).unwrap();
                        }
                    } else {
                        missed_read = missed_read + 1;
                    }
                } else {
                    // log for unannotated splicing junction
                }
            }
        }
    }
    return missed_read;
}