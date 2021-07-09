use std::env;

use coitrees::{COITree, IntervalNode, SortedQuerent};
use std::collections::HashMap;
use std::collections::LinkedList;

extern crate fnv;
use fnv::FnvHashMap;

mod annotations;
mod bam;

fn main() {
    let ann_file_adr: String = env::args().nth(1).expect("missing src");
    let mut transcripts_map: HashMap<String, i32> = HashMap::new();
    let mut transcripts: LinkedList<String> = LinkedList::new();
    let mut tx_lengths: Vec<i32> = Vec::new();
    let trees = annotations::build_tree(&ann_file_adr, &mut transcripts_map, &mut transcripts, &mut tx_lengths);
    // annotations::test_tree(&ann_file_adr, trees.expect("cannot build the tree!"));
    let bam_file_adr: String = env::args().nth(2).expect("missing src");
    bam::read_bamfile(&bam_file_adr, &transcripts_map, &transcripts, &tx_lengths);
}
