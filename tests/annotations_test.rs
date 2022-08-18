// #![deny(warnings)]
// #[macro_use]

use mudskipper::annotation;

use std::collections::HashMap;

#[test]
fn test_gtf_read() {
    let ann_file_adr = "tests/NC_002333.2.gtf".to_string();
    let reader = annotation::read(&ann_file_adr);
    assert!(reader.is_ok());
}

#[test]
fn test_gff_read() {
    let ann_file_adr = "tests/NC_002333.2.gff".to_string();
    let reader = annotation::read(&ann_file_adr);
    assert!(reader.is_ok());
}

#[test]
pub fn test_tree() {
    let ann_file_adr = "tests/NC_002333.2.gff".to_string();
    let mut transcripts_map: HashMap<String, i32> = HashMap::new();
    let mut transcripts: Vec<String> = Vec::new();
    let mut tx_lengths: Vec<i32> = Vec::new();

    let trees = annotation::build_tree(&ann_file_adr, &mut transcripts_map, &mut transcripts, &mut tx_lengths, None).expect("cannot build the tree!");

    let query_file_adr = "tests/NC_002333.2.gff".to_string();
    let reader = annotation::read(&query_file_adr);
    assert!(reader.is_ok());
    let mut failed_count = 0;
    for record in reader.expect("Error reading file.").records() {
        let rec = record.ok().expect("Error reading record.");
        let seqname = rec.seqname().to_string();
        if rec.feature_type() == "exon" {
            let exon_start = *rec.start() as i32;
            let exon_end = *rec.end() as i32;
            if let Some(seqname_tree) = trees.get(&seqname) {
                let countcov = seqname_tree.coverage(exon_start, exon_end);
                let count = countcov.0;
                let cov = countcov.1;
                if count == 0 || cov == 0 {
                    failed_count = failed_count + 1;
                }
            }
        }
    }
    assert_eq!(failed_count, 0, "{} exons are not found in the tree!", failed_count);
}
