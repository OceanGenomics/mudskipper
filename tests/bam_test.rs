// #![deny(warnings)]
// #[macro_use]

extern crate mudskipper;
use mudskipper::annotations;
use mudskipper::bam;

use std::collections::HashMap;

#[test]
pub fn test_read_bamfile() {
    let ann_file_adr = "tests/NC_002333.2.gtf".to_string();
    let mut transcripts_map: HashMap<String, i32> = HashMap::new();
    let mut transcripts: Vec<String> = Vec::new();
    let mut txp_lengths: Vec<i32> = Vec::new();
    let trees = annotations::build_tree(&ann_file_adr, 
                                        &mut transcripts_map,
                                        &mut transcripts,
                                        &mut txp_lengths).expect("cannot build the tree!");
    
    // let bam_file_in = "tests/NC_007112.7.sam".to_string();
    // let bam_file_out = "tests/NC_007112.7.converted.sam".to_string();
    let bam_file_in = "tests/NC_002333.2.sam".to_string();
    let bam_file_out = "tests/NC_002333.2_toTranscriptome.bam".to_string();    
    let missed_count = bam::read_bamfile(&bam_file_in, &bam_file_out, &transcripts, &txp_lengths, &trees);
    let number_of_missed_records = 33;

    // let bam_file_truth = "tests/NC_002333.2_toTranscriptome_truth.bam".to_string();    
    // assert!(compare_bam_files(&bam_file_truth, &bam_file_out), "Some bam records are not converted appropriately.");
    assert_eq!(missed_count, number_of_missed_records, 
                "no transcript for {} alignment records! the correct number is {}.", 
                missed_count,
                number_of_missed_records);
}