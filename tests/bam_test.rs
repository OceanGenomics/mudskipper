// #![deny(warnings)]
// #[macro_use]

use std::str;

use mudskipper::annotation;
use mudskipper::bam;

use std::collections::HashMap;

use rust_htslib::bam::{Read, Reader};


pub fn read_and_process(ann_file_adr: &String,
                        bam_file_in: &String,
                        bam_file_out: &String) -> i32 {
    let mut transcripts_map: HashMap<String, i32> = HashMap::new();
    let mut transcripts: Vec<String> = Vec::new();
    let mut txp_lengths: Vec<i32> = Vec::new();
    let trees = annotation::build_tree(ann_file_adr, 
                                        &mut transcripts_map,
                                        &mut transcripts,
                                        &mut txp_lengths, None).expect("cannot build the tree!");
    let threads = 0;
    let v : Vec<&str> = vec![];
    bam::bam2bam(bam_file_in, bam_file_out, &transcripts, &txp_lengths, &trees, &threads, &200, &v)
}

/*
#[test]
pub fn test_bam2bam() {
    let ann_file_adr = "tests/NC_002333.2.gtf".to_string();
    // let bam_file_in = "tests/NC_007112.7.sam".to_string();
    // let bam_file_out = "tests/NC_007112.7.converted.sam".to_string();
    let bam_file_in = "tests/NC_002333.2.sam".to_string();
    let bam_file_out = "tests/NC_002333.2_toTranscriptome.bam".to_string();    
 
    let number_of_missed_records = 33;
    let missed_count = read_and_process(&ann_file_adr, &bam_file_in, &bam_file_out);

    // let bam_file_truth = "tests/NC_002333.2_toTranscriptome_truth.bam".to_string();    
    // assert!(compare_bam_files(&bam_file_truth, &bam_file_out), "Some bam records are not converted appropriately.");
    assert_eq!(missed_count, number_of_missed_records, 
                "no transcript for {} alignment records! the correct number is {}.", 
                missed_count,
                number_of_missed_records);
}
*/

#[test]
pub fn test_reverse_coordinates_paired() {
    let ann_file_adr = "tests/NC_007112.7_1631.gtf".to_string();
    let bam_file_in = "tests/NC_007112.7_1631.sam".to_string();
    let bam_file_out = "tests/NC_007112.7_1631_toTranscriptome.bam".to_string();

    let missed_count = read_and_process(&ann_file_adr, &bam_file_in, &bam_file_out);
    assert_eq!(missed_count, 0, "no transcripts found!");


    let mut output_bam = Reader::from_path(bam_file_out).unwrap();
    let bam_records = output_bam.records();

    let forward_pos = 321;
    let reverse_pos = 665;

    let mut first = true;
    for rec in bam_records {
        let record = rec.unwrap();
        if first {
            assert_eq!(record.is_paired(), true,  "The alignment record should have been paired!");
            assert_eq!(record.pos(), forward_pos, "Forward pos is wrong! {} {}", record.pos(), forward_pos);
            first = false;
        } else {
            assert_eq!(record.pos(), reverse_pos, "Reverse pos is wrong! {} {}", record.pos(), reverse_pos);
        }

    }
}


#[test]
pub fn test_reverse_coordinates_single() {
    let ann_file_adr = "tests/NC_007112.7_1525.gtf".to_string();
    let bam_file_in = "tests/NC_007112.7_1525.sam".to_string();
    let bam_file_out = "tests/NC_007112.7_1525_toTranscriptome.bam".to_string();

    let missed_count = read_and_process(&ann_file_adr, &bam_file_in, &bam_file_out);
    assert_eq!(missed_count, 0, "no transcripts found!");

    let mut output_bam = Reader::from_path(&bam_file_out).unwrap();
    let bam_records = output_bam.records();
    let output_header = Reader::from_path(&bam_file_out).unwrap();
    let header_view = output_header.header();

    let pos1 = 1401;
    let pos2 = 1111;

    for rec in bam_records {
        let record = rec.unwrap();
        println!("there we are: pos: {}", record.pos());
        let tname_vec = header_view.tid2name(record.tid() as u32).to_vec();
        let tname = str::from_utf8(&tname_vec).unwrap();
        match tname {
            "XM_005159965.4" => assert_eq!(record.pos(), pos1, "Pos1 is wrong! {} {}", record.pos(), pos1),
            "NM_131356.1" => assert_eq!(record.pos(), pos2, "Pos2 is wrong! {} {}", record.pos(), pos2),
            _ => assert!(true, "Target not found!")
        }
    }
}

#[test]
pub fn test_insert_size1() {
    let ann_file_adr = "tests/NC_007120.7_12427.gtf".to_string();
    let bam_file_in = "tests/NC_007120.7_12427.sam".to_string();
    let bam_file_out = "tests/NC_007120.7_12427_toTranscriptome.bam".to_string();

    let missed_count = read_and_process(&ann_file_adr, &bam_file_in, &bam_file_out);
    assert_eq!(missed_count, 0, "no transcripts found!");

    let mut output_bam = Reader::from_path(&bam_file_out).unwrap();
    let bam_records = output_bam.records();

    let insert_size = 389;

    let mut first = true;
    for rec in bam_records {
        let record = rec.unwrap();
        if first {
            assert_eq!(record.is_paired(), true,  "The alignment record should have been paired!");
            assert_eq!(record.insert_size(), insert_size, 
                        "Insert size is wrong! {} {}", record.insert_size(), insert_size);
            first = false;
        } else {
            assert_eq!(record.insert_size(), -insert_size, 
                        "Insert size is wrong! {} {}", record.insert_size(), -insert_size);
        }
    }
}

// different fw and rv length
#[test]
pub fn test_insert_size2() {
    let ann_file_adr = "tests/NC_007114.7_12427.gtf".to_string();
    let bam_file_in = "tests/NC_007114.7_12427.sam".to_string();
    let bam_file_out = "tests/NC_007114.7_12427_toTranscriptome.bam".to_string();

    let missed_count = read_and_process(&ann_file_adr, &bam_file_in, &bam_file_out);
    assert_eq!(missed_count, 0, "no transcripts found!");

    let mut output_bam = Reader::from_path(&bam_file_out).unwrap();
    let bam_records = output_bam.records();

    let insert_size = 264;

    let mut first = true;
    for rec in bam_records {
        let record = rec.unwrap();
        if first {
            assert_eq!(record.is_paired(), true,  "The alignment record should have been paired!");
            assert_eq!(record.insert_size(), insert_size, 
                        "Insert size is wrong! {} {}", record.insert_size(), insert_size);
            first = false;
        } else {
            assert_eq!(record.insert_size(), -insert_size, 
                        "Insert size is wrong! {} {}", record.insert_size(), -insert_size);
        }
    }
}

// different fw and rc length and map to a rc transcript
#[test]
pub fn test_insert_size3() {
    let ann_file_adr = "tests/NC_007131.7_6185.gtf".to_string();
    let bam_file_in = "tests/NC_007131.7_6185.sam".to_string();
    let bam_file_out = "tests/NC_007131.7_6185_toTranscriptome.bam".to_string();

    let missed_count = read_and_process(&ann_file_adr, &bam_file_in, &bam_file_out);
    assert_eq!(missed_count, 0, "no transcripts found!");

    let mut output_bam = Reader::from_path(&bam_file_out).unwrap();
    let bam_records = output_bam.records();

    let insert_size = 455;

    let mut first = true;
    for rec in bam_records {
        let record = rec.unwrap();
        if first {
            assert_eq!(record.is_paired(), true,  "The alignment record should have been paired!");
            assert_eq!(record.insert_size(), insert_size, 
                        "Insert size is wrong! {} {}", record.insert_size(), insert_size);
            first = false;
        } else {
            assert_eq!(record.insert_size(), -insert_size, 
                        "Insert size is wrong! {} {}", record.insert_size(), -insert_size);
        }
    }
}
