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
        // println!("{} {}", res.0, res.1);
        if res.0 != 0 || res.1 != 0 {
            tree.query(range.0, range.1, |node| {
               // println!("query for {} {}:{}",range.0, range.1, node.metadata);
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

pub fn find_ranges_paired(read_pos1: &i32,
                   cigar1: &record::CigarStringView, 
                   new_cigar1: &mut record::CigarString, 
                   read_pos2: &i32, 
                   cigar2: &record::CigarStringView,
                   new_cigar2: &mut record::CigarString) -> Vec<(i32, i32)> {
    let mut ranges1 = find_ranges_single(read_pos1, cigar1, new_cigar1);
    let mut ranges2 = find_ranges_single(read_pos2, cigar2, new_cigar2);
    ranges1.append(&mut ranges2);
    return ranges1;
}

pub fn find_ranges_single(read_pos: &i32, 
                          cigar: &record::CigarStringView,
                          new_cigar_view: &mut record::CigarString) -> Vec<(i32, i32)> {
    let mut ranges = Vec::new();
    let mut curr_range: (i32, i32) = (-1, -1);

    // This variable declares whether a new exon is reached 
    // or the next cigar item belongs to the current exon.
    let mut end_range: bool = true;

    // Set the current position to the base right before the beginning of the alignment
    // after discarding the softclipped bases
    let mut curr_pos = *read_pos - 1 + cigar.leading_softclips() as i32;
    let mut new_cigar : Vec::<record::Cigar> = Vec::<record::Cigar>::new();
    for cigar_item in cigar.iter() {
        // println!("{} {}", cigar_item.char(), cigar_item.len());
        let cigar_char = cigar_item.char();
        let cigar_len = cigar_item.len();
        match cigar_char {
            'M' | 'I' | 'D' => {
                if end_range {
                    curr_range = (curr_pos + 1, -1);
                    end_range = false;
                }
                curr_pos = curr_pos + cigar_len as i32;
                curr_range.1 = curr_pos;
                new_cigar.push(*cigar_item);
            }
            // Observing N means that the rest of 
            // the cigar belongs to another exon
            'N' => {
                end_range = true;
                ranges.push(curr_range);
                curr_pos = curr_pos + cigar_len as i32;
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
    *new_cigar_view = record::CigarString(new_cigar);
    return ranges;
}