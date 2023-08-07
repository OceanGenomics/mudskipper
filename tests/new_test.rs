use std::fs;

use mudskipper::position;



#[test]
fn test_valid_input_data() {
    let input_bam_filename = "/Users/tejastemker/PycharmProjects/Ocean_genomics/mudskipper/tests/input.bam".to_string();

    let mut bqr =  bam2rad_bulk::new(&input_bam_filename, Some(1));

    assert!(bqr.bam2rad_bulk_se().is_ok(),"Mudskipper encountered an error or warning with the input file.");
}


#[test]
fn test_invalid_file_format() {

    let input_file = "tests/unsupported_file.txt".to_string();
    let mut bqr = BAMQueryRecordReader::new(input_file, Some(1));
    let result = bqr.get_next_query_records();
    match result {
        Ok(_) => panic!("Mudskipper should have raised an error for an unsupported file format."),
        Err(error) => {
            assert_eq!(
                error.to_string(),
                "Unsupported file format. The provided file format is not valid or supported."
            );
        }
    }
}


#[test]
fn test_smallest_sam_bam_file() {
    let file_path = "/Users/tejastemker/PycharmProjects/Ocean_genomics/mudskipper/tests/not_name_sorted.sam"; 
    let result = handle_file("/Users/tejastemker/PycharmProjects/Ocean_genomics/mudskipper/tests/not_name_sorted.sam");

    // Assert that no errors or warnings are raised
    assert!(result.is_ok(), "Error occurred while handling the file");
}

fn calculate_mapping_rate(aligned_reads: &[AlignedRead]) -> f64 {
    let total_reads = aligned_reads.len();
    let aligned_count = aligned_reads.iter().filter(|read| read.is_aligned()).count();
    let mut lower_bound = 80.0;
    let mut upper_bound = 100.0;
    let mapping_rate = (aligned_count as f64 / total_reads as f64) * 100.0;

    assert!(mapping_rate >= lower_bound && mapping_rate <= upper_bound,
                "Mapping rate is outside the acceptable range");
}


