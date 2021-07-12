use rust_htslib::bam::record;
use coitrees::{COITree}; //, IntervalNode, SortedQuerent};
use std::collections::HashSet;
// use std::error::Error;

// extern crate fnv;
// use fnv::FnvHashMap;

use crate::annotations;
use annotations::ExonNode;

// type GenericError = Box<dyn Error>;

pub fn find_tid(tree: &COITree<ExonNode, u32>, ranges: &Vec<(i32, i32)>) -> Vec<i32> {
    let mut tids: Vec<i32> = Vec::new();
    let mut tids_set: HashSet<i32> = HashSet::new();
    let mut first = true;
    for range in ranges {
        let res = tree.coverage(range.0, range.1);
        let mut curr_tids: HashSet<i32> = HashSet::new();
        if res.0 != 0 || res.1 != 0 {
            tree.query(range.0, range.1, |node| {
               println!("query for {} {}:{}",range.0, range.1, node.metadata);
               curr_tids.insert(node.metadata.tid);
            });
            if first {
                tids_set = curr_tids;
                first = false;
            } else {
                tids_set = &tids_set & &curr_tids;
            }
            // println!("found coverage: {:?}", res)
        }
    }
    for tid in tids_set {
        tids.push(tid);
    }
    return tids;
}

pub fn find_ranges(read_pos: &i32, cigar: record::CigarStringView) -> Vec<(i32, i32)> {
    let mut ranges = Vec::new();
    let mut curr_range: (i32, i32) = (-1, -1);
    let mut end_range: bool = true;
    let mut curr_pos = *read_pos - 1 + cigar.leading_softclips() as i32;
    for cigar_item in cigar.iter() {
        let cigar_char = cigar_item.char();
        let cigar_len = cigar_item.len();
        match cigar_char {
            'M' | 'I' | 'D' => {
                if end_range {
                    curr_range = (curr_pos + 1, -1);
                //    curr_pos = curr_pos + cigar_len as i32;
                //    curr_range.1 = curr_pos;
                }
                curr_pos = curr_pos + cigar_len as i32;
                curr_range.1 = curr_pos;
            }
            'N' => {
                end_range = true;
                ranges.push(curr_range);
                curr_pos = curr_pos + cigar_len as i32;
            }
            /*'I' => { 
                if end_range {
                    curr_range = (curr_pos + 1 , -1);
                }
                curr_pos = curr_pos + cigar_len as i32;
                curr_range.1 = curr_pos;
            }
            'D' => {
                if end_range {
                    curr_range = (curr_pos + 1 , -1);
                }
                curr_pos = curr_pos + cigar_len as i32;
                curr_range.1 = curr_pos;
            }*/
            'S' => { 
                // curr_pos = curr_pos + cigar_len as i32;
            }
            _ => println!("Unexpected cigar char! {}", cigar_char)
        }
    }
    ranges.push(curr_range);
    return ranges;
}