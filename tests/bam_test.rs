// #![deny(warnings)]
// #[macro_use]

extern crate mudskipper;
use mudskipper::annotations;
use mudskipper::bam;
use mudskipper::intersection;

use std::collections::HashMap;

#[test]
pub fn test_read_bamfile() {
    let ann_file_adr = "tests/NC_007112.7.gtf".to_string();
    let mut transcripts_map: HashMap<String, i32> = HashMap::new();
    let mut transcripts: Vec<String> = Vec::new();
    let mut txp_lengths: Vec<i32> = Vec::new();
    let trees = annotations::build_tree(&ann_file_adr, 
                                        &mut transcripts_map,
                                        &mut transcripts,
                                        &mut txp_lengths).expect("cannot build the tree!");
    
    // let bam_file_in = "tests/NC_007112.7.sam".to_string();
    // let bam_file_out = "tests/NC_007112.7.converted.sam".to_string();
    let bam_file_in = "/home/mohsen/d_rerio/star_NC_007112.7/Aligned.out.sam".to_string();
    let bam_file_out = "/home/mohsen/d_rerio/star_NC_007112.7/Aligned.out2.sam".to_string();    
    let missed_count = bam::read_bamfile(&bam_file_in, &bam_file_out, &transcripts, &txp_lengths, &trees);
    assert_eq!(missed_count, 0, "no transcript for {} alignment records!", missed_count);
}