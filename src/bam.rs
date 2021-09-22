use coitrees::{COITree};

use rust_htslib::bam::{Format, Header, Read, Reader, Writer, header, record};

use crate::annotation;
use annotation::ExonNode;

extern crate bio_types;

use crate::convert;

extern crate fnv;
use fnv::FnvHashMap;

use log::{debug, error};

// pub fn bam2bam(input_bam_filename: &String, 
//                output_bam_filename: &String,
//                transcripts: &Vec<String>,
//                txp_lengths: &Vec<i32>,
//                trees: &FnvHashMap::<String, COITree<ExonNode, u32>>,
//                threads_count: &usize,
//                max_softlen: &usize) -> i32 {
//     let mut input_bam = Reader::from_path(input_bam_filename).unwrap();
//     // let bam_records = input_bam.records();

//     // let header_ = Header::from_template(input_bam.header());
//     // let header_text = String::from_utf8(header_.to_bytes().to_owned()).unwrap().;

//     let input_bam_header = Reader::from_path(input_bam_filename).unwrap();
//     let header_view = input_bam_header.header();
    
//     let mut new_header = Header::new();
//     // new_header.push_record(header::HeaderRecord::new(b"HD").push_tag(b"VN", &"1.4"));

//     let mut counter = 0;
//     for tid in transcripts.iter() {
//         let txp_len = txp_lengths[counter];
//         new_header.push_record(header::HeaderRecord::new(b"SQ").push_tag(b"SN", &tid).push_tag(b"LN", &txp_len));
//         counter += 1;
//     }

//     let mut output_bam = Writer::from_path(output_bam_filename, &new_header, Format::Bam).unwrap();

//     if *threads_count >= 2 {
//         let threads_count_half = threads_count / 2;
//         println!("thread count: {}", threads_count_half);
//         input_bam.set_threads(threads_count_half).expect("Failed to set number of BAM reading threads.");
//         output_bam.set_threads(threads_count_half).expect("Failed to set number of BAM writing threads.");
//     }

//     let mut n = 0;
//     let mut missed_read = 0;
//     let mut first_in_pair = true;
//     let mut first_record: record::Record = record::Record::new();
//     let mut new_cigar: record::CigarString = record::CigarString(vec![record::Cigar::Match(100)]);
//     let mut first_new_cigar: record::CigarString = record::CigarString(vec![record::Cigar::Match(100)]);
//     let mut len1 = 0;
//     let mut len2 = 0;
//     for rec in input_bam.records() {
//         n = n + 1;
//         let record = rec.unwrap();
//         let qname = String::from_utf8(record.qname().to_vec()).unwrap();
//         // let mut check = false;
//         // if qname == "J00153:160:HMTNMBBXX:6:1101:14905:6185" {
//         //     check  = true;
//         // }
//         debug!("qname: {}", qname);
//         let mut long_softclip = false;
//         if !record.is_paired() || record.is_mate_unmapped() {
//             let ranges = convert::find_ranges_single(&(record.pos() as i32),
//                                                           &record.cigar(),
//                                                           &mut new_cigar,
//                                                           &mut len1,
//                                                           &mut long_softclip,
//                                                           &max_softlen);
//             let genome_tname = String::from_utf8(header_view.tid2name(record.tid() as u32).to_vec())
//                                                 .expect("cannot find the tname!");
//             if let Some(tree) = trees.get(&genome_tname) {
//                 let tids = convert::find_tid(&tree, &ranges);
//                 if long_softclip {
//                     debug!("The softclip length is too long!");
//                     continue;
//                 }
//                 if tids.len() > 0 {
//                     for (tid, pos_strand) in tids.iter() {
//                         debug!("{} {}", tid, transcripts[*tid as usize]);
//                         debug!("{}", tid);

//                         let mut record_ = record.clone();
//                         let mut pos = 0;
//                         if pos_strand.1 == Strand::Forward {
//                             pos = record.pos() - (pos_strand.0 as i64);
//                         } else if pos_strand.1 == Strand::Reverse { 
//                             pos = (pos_strand.0 as i64) - record.pos() - len1 as i64;
//                             if record.is_reverse() {
//                                 record_.unset_reverse();
//                             } else {
//                                 record_.set_reverse();
//                             }
//                         }

//                         record_.set(record.qname(),
//                                     Some(&new_cigar),
//                                     &record.seq().as_bytes(),
//                                     record.qual());
//                         record_.set_tid(*tid);
//                         record_.set_pos(pos);

//                         output_bam.write(&record_).unwrap();
//                     }
//                 } else {
//                     missed_read = missed_read + 1;
//                 }
//             } else {
//                 // log for unannotated splicing junction
//             }
//             debug!("{} {}", ranges[0].0, ranges[0].1);
//             first_in_pair = true;
//         } else {
//             if first_in_pair {
//                 first_in_pair = false;
//                 first_record = record;
//             } else {
//                 // if check {
//                 //     println!("Here is reached!");
//                 // }
//                 first_in_pair = true;
//                 /* let ranges = convert::find_tids_paired(&(first_record.pos() as i32), 
//                                                               &first_record.cigar(),
//                                                               &mut first_new_cigar,
//                                                               &(record.pos() as i32), 
//                                                               &record.cigar(), 
//                                                               &mut new_cigar); */
//                 let genome_tname = String::from_utf8(header_view.tid2name(record.tid() as u32)
//                                                                 .to_vec()).expect("cannot find the tname!");
//                 if let Some(tree) = trees.get(&genome_tname) {
//                     let tids = convert::find_tids_paired(&tree,
//                                                                 &(first_record.pos() as i32), 
//                                                                 &first_record.cigar(),
//                                                                 &mut first_new_cigar,
//                                                                 &mut len1,
//                                                                 &(record.pos() as i32), 
//                                                                 &record.cigar(), 
//                                                                 &mut new_cigar,
//                                                                 &mut len2,
//                                                                 &mut long_softclip,
//                                                                 &max_softlen);
//                     if long_softclip {
//                         debug!("The softclip length is too long!");
//                         continue;
//                     }
//                     debug!("{}: {}", first_record.cigar(), first_record.cigar().len());
//                     debug!("{}: {}", record.cigar(), record.cigar().len());
//                     // if check {
//                     //     println!("{}: {}", first_record.cigar(), first_record.cigar().len());
//                     //     println!("{}: {}", record.cigar(), record.cigar().len());
//                     // }
//                     if tids.len() > 0 {
//                         for (tid, pos_strand) in tids.iter() {
//                             let mut first_record_ = first_record.clone();
//                             let mut second_record_ = record.clone();

//                             let mut first_pos = 0;
//                             let mut second_pos = 0;
//                             if pos_strand.0.1 == Strand::Forward {
//                                 first_pos = first_record.pos() - (pos_strand.0.0 as i64);
//                                 debug!("first_pos:{} - pos:{} = {}",
//                                         first_record.pos(), pos_strand.0.0, first_pos);
//                                 second_pos = record.pos() - (pos_strand.1.0 as i64); 
//                                 debug!("second_pos:{} - pos:{} = {}",
//                                     record.pos(), pos_strand.1.0, second_pos);
//                             } else if pos_strand.0.1 == Strand::Reverse{ 
//                                 first_pos = (pos_strand.0.0 as i64) - first_record.pos() - len1 as i64;
//                                 debug!("pos:{} - first_pos:{} - len1:{} = {}",
//                                         pos_strand.0.0, first_record.pos(), len1, first_pos);
//                                 second_pos = (pos_strand.1.0 as i64) - record.pos() - len2 as i64;
//                                 debug!("pos:{} - second_pos:{} - len2:{} = {}",
//                                         pos_strand.1.0, record.pos(), len2, second_pos);
//                                 if first_record.is_reverse() {
//                                     debug!("first_record is reversed");
//                                     first_record_.unset_reverse();
//                                     first_record_.set_mate_reverse();
//                                 } else {
//                                     debug!("first_record is not reversed");
//                                     first_record_.set_reverse();
//                                     first_record_.unset_mate_reverse();
//                                 }
//                                 if record.is_reverse() {
//                                     debug!("second_record is reversed");
//                                     second_record_.unset_reverse();
//                                     second_record_.set_mate_reverse();
//                                 } else {
//                                     debug!("second_record is not reversed");
//                                     second_record_.set_reverse();
//                                     second_record_.unset_mate_reverse();
//                                 }

//                             }
//                             let first_read_len: i64 = first_record.seq().len() as i64;
//                             let second_read_len: i64 = record.seq().len() as i64;
//                             let first_length: i64;
//                             let second_length: i64;
//                             if pos_strand.0.1 == Strand::Forward {
//                                 first_length = if !first_record.is_reverse() {
//                                                     second_pos - first_pos + second_read_len } else {
//                                                     second_pos - first_pos + first_read_len };
//                                 second_length = if !record.is_reverse() {
//                                                     first_pos - second_pos - first_read_len } else {
//                                                     first_pos - second_pos - second_read_len };
//                             }
//                             else {
//                                 first_length = if !first_record.is_reverse() {
//                                     second_pos - first_pos - first_read_len } else {
//                                     second_pos - first_pos - second_read_len };
//                                 second_length = if !record.is_reverse() {
//                                     first_pos - second_pos + second_read_len } else {
//                                     first_pos - second_pos + first_read_len };
//                             }
//                             debug!("{} {}", tid, transcripts[*tid as usize]);
//                             debug!("first_pos:{} second_pos:{} len1:{} len2:{} first_length:{} second_length:{}",
//                                     first_pos, second_pos, first_read_len, second_read_len, first_length, second_length);
//                             debug!("first_record.is_reverse():{}", first_record.is_reverse());
//                             debug!("record.is_reverse():{}", record.is_reverse());
//                             // if check {
//                             //     println!("{} {}", tid, transcripts[*tid as usize]);
//                             //     println!("first_pos:{} second_pos:{} len1:{} len2:{} first_length:{} second_length:{}",
//                             //             first_pos, second_pos, first_read_len, second_read_len, first_length, second_length);
//                             //     println!("first_record.is_reverse():{}", first_record.is_reverse());
//                             //     println!("record.is_reverse():{}", record.is_reverse());    
//                             // }
//                             first_record_.set(first_record.qname(), 
//                                               Some(&first_new_cigar), 
//                                               &first_record.seq().as_bytes(), 
//                                               first_record.qual());
//                             first_record_.set_tid(*tid);
//                             first_record_.set_mtid(*tid);

//                             first_record_.set_pos(first_pos);
//                             first_record_.set_mpos(second_pos);
//                             first_record_.set_insert_size(first_length);
//                             // first_record_.remove_aux("AS".as_bytes());

//                             second_record_.set(record.qname(), Some(&new_cigar), 
//                                                &record.seq().as_bytes(), 
//                                                record.qual());
//                             second_record_.set_tid(*tid);
//                             second_record_.set_mtid(*tid);

//                             second_record_.set_pos(second_pos);
//                             second_record_.set_mpos(first_pos);
//                             second_record_.set_insert_size(second_length);
//                             // second_record_.remove_aux("AS".as_bytes());

//                             if first_pos > second_pos {
//                                 if first_pos + first_read_len > txp_lengths[*tid as usize].into() || second_pos < 0 {
//                                     continue;
//                                 }
//                                 output_bam.write(&second_record_).unwrap();
//                                 output_bam.write(&first_record_).unwrap();    
//                             } else {
//                                 if second_pos + second_read_len > txp_lengths[*tid as usize].into() || first_pos < 0 {
//                                     continue;
//                                 }
//                                 output_bam.write(&first_record_).unwrap();
//                                 output_bam.write(&second_record_).unwrap();
//                             }
//                         }
//                     } else {
//                         missed_read = missed_read + 1;
//                         debug!("missed read!");
//                     }
//                 } else {
//                     // log for unannotated splicing junction
//                 }
//             }
//         }
        
//         // if n == 4 {
//         //     break;
//         // }
//     }
//     // output_bam.flush();
//     return missed_read;
// }

pub fn bam2bam(input_bam_filename: &String, 
                   output_bam_filename: &String,
                   transcripts: &Vec<String>,
                   txp_lengths: &Vec<i32>,
                   trees: &FnvHashMap::<String, COITree<ExonNode, u32>>,
                   threads_count: &usize,
                   max_softlen: &usize) {
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
    let mut mate_wanted = false;
    let mut first_record: record::Record = record::Record::new();
    for rec in input_bam.records() {
        n = n + 1;
        let record = rec.unwrap();
        let qname = String::from_utf8(record.qname().to_vec()).unwrap();

        debug!("qname: {}", qname);
        if !record.is_paired() || record.is_mate_unmapped() {
            let txp_records = convert::convert_single_end(&record, 
                                           header_view,
                                           transcripts,
                                           trees,
                                           max_softlen);
            for txp_rec in txp_records.iter() {
                output_bam.write(txp_rec).unwrap();
            }
            mate_wanted = false;
        } else {
            if mate_wanted {
                // this is the second read in pair... PROCESS!
                let txp_records = convert::convert_paired_end(&first_record,
                                                                        &record,
                                                                        header_view,
                                                                        transcripts,
                                                                        txp_lengths,
                                                                        trees,
                                                                        max_softlen);
                for txp_rec in txp_records.iter() {
                    output_bam.write(txp_rec).unwrap();
                }
                mate_wanted = false;
            } else {
                // this is the first read in pair... save and wait for the second read in the pair
                mate_wanted = true;
                first_record = record;
            }
        }
    }
}

pub fn bam2bam_tags(input_bam_filename: &String, 
                   output_bam_filename: &String,
                   transcripts: &Vec<String>,
                   txp_lengths: &Vec<i32>,
                   trees: &FnvHashMap::<String, COITree<ExonNode, u32>>,
                   threads_count: &usize,
                   max_softlen: &usize, 
                   required_tags: &Vec<&str>) {
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
    let mut mate_wanted = false;
    let mut first_record: record::Record = record::Record::new();
    for rec in input_bam.records() {
        n = n + 1;
        let record = rec.unwrap();
        let qname = String::from_utf8(record.qname().to_vec()).unwrap();
        for tag in required_tags.iter() {
            if record.aux(tag.as_bytes()).is_err() {
                error!("Could not find {} tag for read {}", tag, qname);
                panic!("Some required tags do not exist!");
            }
        }

        debug!("qname: {}", qname);
        if !record.is_paired() || record.is_mate_unmapped() {
            let txp_records = convert::convert_single_end(&record, 
                                           header_view,
                                           transcripts,
                                           trees,
                                           max_softlen);
            for txp_rec in txp_records.iter() {
                output_bam.write(txp_rec).unwrap();
            }
            mate_wanted = false;
        } else {
            if mate_wanted {
                // this is the second read in pair... PROCESS!
                let txp_records = convert::convert_paired_end(&first_record,
                                                                        &record,
                                                                        header_view,
                                                                        transcripts,
                                                                        txp_lengths,
                                                                        trees,
                                                                        max_softlen);
                for txp_rec in txp_records.iter() {
                    output_bam.write(txp_rec).unwrap();
                }
                mate_wanted = false;
            } else {
                // this is the first read in pair... save and wait for the second read in the pair
                mate_wanted = true;
                first_record = record;
            }
        }
    }
}