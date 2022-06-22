use mudskipper::query_bam_records::BAMQueryRecordReader;
use mudskipper::position;
use std::panic;

#[test]
#[should_panic]
fn test_missing_mate() {
    // setup the input BAM file
    let input_bam_filename = "tests/missing_mate.sam".to_string();
    // setup the BAM query record reader
    let mut bqr = BAMQueryRecordReader::new(&input_bam_filename, Some(1));
    bqr.get_next_query_records();
}

#[test]
#[should_panic]
fn test_missing_mate2() {
    // setup the input BAM file
    let input_bam_filename = "tests/missing_mate2.sam".to_string();
    // setup the BAM query record reader
    let mut bqr = BAMQueryRecordReader::new(&input_bam_filename, Some(1));
    bqr.get_next_query_records();
}

#[test]
//#[should_panic]
fn test_name_order() {
    // setup the input BAM file
    let input_bam_filename = "tests/not_name_sorted.sam".to_string();
    // setup the BAM query record reader 
    let result = panic::catch_unwind(||{
        let mut bqr = BAMQueryRecordReader::new(&input_bam_filename, Some(1));
        let _next_query = bqr.get_next_query_records();
    });
    assert!(result.is_err());

    // run shuffle
    // bam file is name sorted after shuffle
    position::depositionify_bam(&input_bam_filename, &"tests/name_sorted.bam", 8 * 1024 * 1024, 1);
    let sorted_bam_filename = "tests/name_sorted.bam".to_string();
    let result = panic::catch_unwind(||{
        let mut bqr = BAMQueryRecordReader::new(&sorted_bam_filename, Some(1));
        let _next_query = bqr.get_next_query_records();
    });
    assert!(result.is_ok());
    
}
