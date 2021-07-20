use coitrees::{COITree}; //, IntervalNode, SortedQuerent};

use rust_htslib::bam::{Format, Header, Read, Reader, Writer, header, record}; //, HeaderView};
// use std::collections::HashMap;
// use std::collections::LinkedList;

use crate::annotations;
use annotations::ExonNode;

use crate::intersection;

extern crate fnv;
use fnv::FnvHashMap;

pub fn read_bamfile(input_bam_filename: &String, 
                    output_bam_filename: &String,
                    transcripts: &Vec<String>,
                    txp_lengths: &Vec<i32>,
                    trees: &FnvHashMap::<String, COITree<ExonNode, u32>>) -> i32 {
    let mut input_bam = Reader::from_path(input_bam_filename).unwrap();
    let bam_records = input_bam.records();

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

    // input_bam.set_threads(2).expect("Failed to set number of BAM reading threads to 2.");
    // output_bam.set_threads(5).expect("Failed to set number of BAM writing threads to 4.");

    let mut n = 0;
    let mut missed_read = 0;
    let mut first_in_pair = true;
    let mut first_record : record::Record = record::Record::new();
    let mut new_cigar : record::CigarString = record::CigarString(vec![record::Cigar::Match(100)]);
    let mut first_new_cigar : record::CigarString = record::CigarString(vec![record::Cigar::Match(100)]);
        for rec in bam_records {
        n = n + 1;
        let record = rec.unwrap();
        if !record.is_paired() {
            let ranges = intersection::find_ranges_single(&(record.pos() as i32), &record.cigar(), &mut new_cigar);
            let genome_tname = String::from_utf8(header_view.tid2name(record.tid() as u32).to_vec()).expect("cannot find the tname!");
            if let Some(tree) = trees.get(&genome_tname) {
                let tids = intersection::find_tid(&tree, &ranges);
                // println!("here is reached. {}", tids.len());
                if tids.len() > 0 {
                    // println!("\n {} {} {:?}\n", record.pos(), record.cigar());
                    for (tid, pos) in tids.iter() {
                        // println!("{} {}", tid, transcripts[*tid as usize]);
                        // println!("{}", tid);

                        let mut record_ = record.clone(); //Record::new();
                        record_.set(record.qname(), Some(&new_cigar), &record.seq().as_bytes(), record.qual());
                        record_.set_tid(*tid);
                        record_.set_pos(record.pos() - (*pos as i64));

                        output_bam.write(&record_).unwrap();   
                    }
                    // println!("done!");
                } else {
                    // println!("\n {} {} {:?}\n", record.pos(), record.cigar());
                    missed_read = missed_read + 1;
                }
            } else {
                // log for unannotated splicing junction
            }
            // println!("{} {}", ranges[0].0, ranges[0].1);
            
        } else {
            if first_in_pair {
                first_in_pair = false;
                first_record = record;
            } else {
                first_in_pair = true;
                let ranges = intersection::find_ranges_paired(&(first_record.pos() as i32), 
                                                              &first_record.cigar(),
                                                              &mut first_new_cigar,
                                                              &(record.pos() as i32), 
                                                              &record.cigar(), 
                                                              &mut new_cigar);
                let genome_tname = String::from_utf8(header_view.tid2name(record.tid() as u32).to_vec()).expect("cannot find the tname!");
                if let Some(tree) = trees.get(&genome_tname) {
                    let tids = intersection::find_tid(&tree, &ranges);
                    // println!("here is reached. {}", tids.len());
                    if tids.len() > 0 {
                        // println!("\n {} {} {:?}\n", record.pos(), record.cigar());
                        for (tid, pos) in tids.iter() {
                            // println!("{} {}", tid, transcripts[*tid as usize]);
                            // println!("{}", tid);

                            let mut first_record_ = first_record.clone(); //Record::new();
                            first_record_.set(first_record.qname(), Some(&first_new_cigar), &first_record.seq().as_bytes(), first_record.qual());
                            first_record_.set_tid(*tid);
                            first_record_.set_pos(first_record.pos() - (*pos as i64));
                            // first_record_.set_pos();
                            // first_record_.set_mpos();
                            
                            output_bam.write(&first_record_).unwrap();

                            let mut second_record_ = record.clone(); //Record::new();
                            second_record_.set(record.qname(), Some(&new_cigar), &record.seq().as_bytes(), record.qual());
                            second_record_.set_tid(*tid);
                            // second_record_.set_pos();
                            // second_record_.set_mpos();

                            output_bam.write(&second_record_).unwrap();                            
                        }
                        // println!("done!");
                    } else {
                        // println!("\n {} {} {:?}\n", record.pos(), record.cigar());
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