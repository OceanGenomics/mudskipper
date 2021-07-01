use std::env;

use coitrees::{COITree, IntervalNode, SortedQuerent};

extern crate fnv;
use fnv::FnvHashMap;

mod annotations;

fn main() {
    let ann_file_adr: String = env::args().nth(1).expect("missing src");
    let trees = annotations::build_tree(&ann_file_adr);
    annotations::test_tree(&ann_file_adr, trees.expect("cannot build the tree!"));
}
