use std::env;

// use coitrees::{COITree, IntervalNode, SortedQuerent};
use std::collections::HashMap;
// use std::collections::LinkedList;

// extern crate fnv;
// use fnv::FnvHashMap;

mod annotations;
mod intersection;
mod bam;

fn main() {
    let ann_file_adr: String = env::args().nth(1).expect("missing src");
    let mut transcripts_map: HashMap<String, i32> = HashMap::new();
    let mut transcripts: Vec<String> = Vec::new();
    let mut txp_lengths: Vec<i32> = Vec::new();
    let trees = annotations::build_tree(&ann_file_adr, 
                                        &mut transcripts_map,
                                        &mut transcripts,
                                        &mut txp_lengths).expect("cannot build the tree!");
    let bam_file_in: String = env::args().nth(2).expect("missing src");
    let bam_file_out: String = env::args().nth(3).expect("missing src");
    bam::read_bamfile(&bam_file_in, &bam_file_out, &transcripts, &txp_lengths, &trees);
}