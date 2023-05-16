use std::fs;

use mudskipper::position;
use mudskipper::query_bam_records::BAMQueryRecordReader;

#[test]
fn test_missing_mate() {
    // setup the input BAM file
    let input_bam_filename1 = "tests/missing_mate.sam".to_string();
    let input_bam_filename2 = "tests/missing_mate2.sam".to_string();
    let input_bam_filename3 = "tests/missing_mate3.sam".to_string();
    // setup the BAM query record reader
    let mut bqr1 = BAMQueryRecordReader::new(&input_bam_filename1, Some(1));
    assert_eq!(bqr1.get_next_query_records().is_err(), true);
    let mut bqr2 = BAMQueryRecordReader::new(&input_bam_filename2, Some(1));
    assert_eq!(bqr2.get_next_query_records().is_err(), true);
    let mut bqr3 = BAMQueryRecordReader::new(&input_bam_filename3, Some(1));
    assert_eq!(bqr3.get_next_query_records().is_err(), false);
}

#[test]
fn test_name_order() {
    // setup the input BAM file
    let input_bam_filename = "tests/not_name_sorted.sam".to_string();
    // setup the BAM query record reader
    let mut bqr = BAMQueryRecordReader::new(&input_bam_filename, Some(1));
    assert_eq!(bqr.get_next_query_records().is_err(), true);

    // run shuffle
    // bam file is name sorted after shuffle
    position::depositionify_bam(&input_bam_filename, &"tests/name_sorted.bam", 8 * 1024 * 1024, 1);
    let sorted_bam_filename = "tests/name_sorted.bam".to_string();
    let mut bqr = BAMQueryRecordReader::new(&sorted_bam_filename, Some(1));
    assert_eq!(bqr.get_next_query_records().is_ok(), true);
    // remove the name sorted file
    fs::remove_file(sorted_bam_filename).unwrap();
}
