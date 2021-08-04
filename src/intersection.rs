use rust_htslib::bam::record;
use coitrees::{COITree}; //, IntervalNode, SortedQuerent};
use std::collections::{HashSet, HashMap};
// use std::error::Error;

extern crate bio_types;
use bio_types::strand::Strand;

use crate::annotations;
use annotations::ExonNode;

// type GenericError = Box<dyn Error>;

pub fn find_tid(tree: &COITree<ExonNode, u32>, ranges: &Vec<(i32, i32)>) ->HashMap<i32, (i32, Strand)> {
    let mut tids: HashMap<i32, (i32, Strand)> = HashMap::new();
    let mut tids_set: HashSet<i32> = HashSet::new();
    let mut tid_pos: HashMap<i32, (i32, Strand)> = HashMap::new();
    let mut first = true;
    let mut last = false;
    let mut counter = 1;
    for range in ranges {
        if counter == ranges.len() {
            last = true;
        }
        counter = counter + 1;
        let res = tree.coverage(range.0, range.1);
        let mut curr_tids: HashSet<i32> = HashSet::new();
        println!("range: {} {}", res.0, res.1);
        println!("first: {}, last: {}", first, last);
        if res.0 != 0 || res.1 != 0 {
            tree.query(range.0, range.1, |node| {
                println!("query for {} {}:{}",range.0, range.1, node.metadata);
                /* println!("start: {} - tpos_start: {} = {}",
                        node.metadata.end, node.metadata.tpos_start, 
                        node.metadata.end - node.metadata.tpos_start); */

                curr_tids.insert(node.metadata.tid);
                if first && node.metadata.strand == Strand::Forward {
                    println!("inserting: {}", node.metadata.tid);
                    tid_pos.insert(node.metadata.tid, 
                        (node.metadata.start - node.metadata.tpos_start - 1, Strand::Forward));
                } else if last && node.metadata.strand == Strand::Reverse {
                    println!("inserting: {}", node.metadata.tid);
                    tid_pos.insert(node.metadata.tid, 
                        (node.metadata.end + node.metadata.tpos_start, Strand::Reverse));
                }
            });
            // println!("found coverage: {:?}", res)
        }
        if first {
            tids_set = curr_tids;
            first = false;
        } else {
            tids_set = &tids_set & &curr_tids;
        }
    }
    for tid in tids_set {
        println!("querying: {}", tid);
        tids.insert(tid, tid_pos[&tid]);
    }
    return tids;
}

pub fn find_tids_paired(tree: &COITree<ExonNode, u32>,
                        read_pos1: &i32,
                        cigar1: &record::CigarStringView, 
                        new_cigar1: &mut record::CigarString, 
                        len1: &mut i32,
                        read_pos2: &i32, 
                        cigar2: &record::CigarStringView,
                        new_cigar2: &mut record::CigarString,
                        len2: &mut i32) -> HashMap<i32, ((i32, Strand), (i32, Strand))> {
    let mut tid_pos: HashMap<i32, ((i32, Strand), (i32, Strand))> = HashMap::new();
    println!("read1: {} {}", read_pos1, cigar1);
    let ranges1 = find_ranges_single(read_pos1, cigar1, new_cigar1, len1);   
    let tids1 = find_tid(&tree, &ranges1);
    println!("read2: {} {}", read_pos2, cigar2);
    let ranges2 = find_ranges_single(read_pos2, cigar2, new_cigar2, len2);
    let tids2 = find_tid(&tree, &ranges2);
    for (tid, pos1) in tids1.iter() {
        if tids2.contains_key(tid) {
            tid_pos.insert(*tid, (*pos1, tids2[tid]));
        }
    }
    return tid_pos;
}

pub fn find_ranges_single(read_pos: &i32, 
                          cigar: &record::CigarStringView,
                          new_cigar_view: &mut record::CigarString,
                          len: &mut i32) -> Vec<(i32, i32)> {
    let mut ranges = Vec::new();
    let mut curr_range: (i32, i32) = (-1, -1);

    // This variable declares whether a new exon is reached 
    // or the next cigar item belongs to the current exon.
    let mut end_range: bool = true;

    // Set the current position to the base right before the beginning of the alignment
    // after discarding the softclipped bases
    let mut curr_pos = *read_pos - 1 + (cigar.leading_softclips() as i32);
    let mut new_cigar : Vec::<record::Cigar> = Vec::<record::Cigar>::new();
    *len = 0;
    for cigar_item in cigar.iter() {
        println!("cigar item: {} {}", cigar_item.char(), cigar_item.len());
        let cigar_char = cigar_item.char();
        let mut cigar_len = cigar_item.len();
        match cigar_char {
            'M' | 'I' | 'D' => {
                if new_cigar.len() == 0 || new_cigar.last().unwrap().char() != cigar_char {
                    new_cigar.push(*cigar_item);
                } else {
                    cigar_len = cigar_len + new_cigar.last().unwrap().len();
                    new_cigar.pop();
                    match cigar_char {
                        'M' => new_cigar.push(record::Cigar::Match(cigar_len)),
                        'I' => new_cigar.push(record::Cigar::Ins(cigar_len)),
                        'D' => new_cigar.push(record::Cigar::Del(cigar_len)),
                        _ => println!("Warning! unexpected cigar item: {}", cigar_char)
                    }
                }                
                if cigar_char != 'I' {
                    *len = *len + cigar_len as i32; 
                }

                if end_range {
                    curr_range = (curr_pos + 1, -1);
                    end_range = false;
                }
                curr_pos = curr_pos + cigar_len as i32;
                curr_range.1 = curr_pos;
            }
            // Observing N means that the rest of 
            // the cigar belongs to another exon
            'N' => {
                end_range = true;
                ranges.push(curr_range);
                println!("pushing {} {}", curr_range.0, curr_range.1);
                curr_pos = curr_pos + cigar_len as i32;
                *len = *len + cigar_len as i32;
            }
            // Soft clipped bases at the beginning are
            // taken care of before beginning the loop
            // and the ending softclipped bases should
            // be ignored.
            'S' => {
                new_cigar.push(*cigar_item);
            }
            _ => println!("Unexpected cigar char! {}", cigar_char)
        }
        // println!("{} {}", curr_range.0, curr_range.1);
    }
    ranges.push(curr_range);
    println!("pushing {} {}", curr_range.0, curr_range.1);
    *new_cigar_view = record::CigarString(new_cigar);
    return ranges;
}