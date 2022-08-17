use std::collections::HashMap;

use bio_types::strand::Strand;
use mudskipper::{convert::{convert_cigar_overhang, TranscriptInfo}, annotation, bam};
use rust_htslib::bam::record::{CigarString, Cigar};
use rust_htslib::bam::{Read, Reader};

#[test]
fn test_overhang_cigar() {
    // 100M
    let cigar1 = CigarString(vec![Cigar::Match(100)]);
    // 70M2D28M
    let cigar2 = CigarString(vec![Cigar::Match(70), Cigar::Del(2), Cigar::Match(28)]);
    // 9M1D80M1D9M
    let cigar3 = CigarString(vec![Cigar::Match(9), Cigar::Del(1), Cigar::Match(80), Cigar::Del(1), Cigar::Match(9)]);
    // 7M1D93M
    let cigar4 = CigarString(vec![Cigar::Match(7), Cigar::Del(1), Cigar::Match(93)]);
    // 93M1D7M
    let cigar5 = CigarString(vec![Cigar::Match(93), Cigar::Del(1), Cigar::Match(7)]);
    // 7M1I92M
    let cigar6 = CigarString(vec![Cigar::Match(7), Cigar::Ins(1), Cigar::Match(92)]);
    // 92M1I7M
    let cigar7 = CigarString(vec![Cigar::Match(92), Cigar::Ins(1), Cigar::Match(7)]);
    // 3S97M
    let cigar8 = CigarString(vec![Cigar::SoftClip(3), Cigar::Match(97)]);
    // 97M3S
    let cigar9 = CigarString(vec![Cigar::Match(97), Cigar::SoftClip(3)]);
    
    let tinfo = TranscriptInfo {
        pos: 1,
        strand: Strand::Forward,
        left_overhang: 10,
        right_overhang: 10,
    };

    let converted_cigar1 = convert_cigar_overhang(&cigar1, &tinfo);
    let converted_cigar2 = convert_cigar_overhang(&cigar2, &tinfo);
    let converted_cigar3 = convert_cigar_overhang(&cigar3, &tinfo);
    let converted_cigar4 = convert_cigar_overhang(&cigar4, &tinfo);
    let converted_cigar5 = convert_cigar_overhang(&cigar5, &tinfo);
    let converted_cigar6 = convert_cigar_overhang(&cigar6, &tinfo);
    let converted_cigar7 = convert_cigar_overhang(&cigar7, &tinfo);
    let converted_cigar8 = convert_cigar_overhang(&cigar8, &tinfo);
    let converted_cigar9 = convert_cigar_overhang(&cigar9, &tinfo);

    assert_eq!(format!("{}", converted_cigar1), "10S80M10S");
    assert_eq!(format!("{}", converted_cigar2), "10S60M2D18M10S");
    assert_eq!(format!("{}", converted_cigar3), "10S78M10S");
    assert_eq!(format!("{}", converted_cigar4), "10S80M10S");
    assert_eq!(format!("{}", converted_cigar5), "10S80M10S");
    assert_eq!(format!("{}", converted_cigar6), "10S80M10S");
    assert_eq!(format!("{}", converted_cigar7), "10S80M10S");
    assert_eq!(format!("{}", converted_cigar8), "10S80M10S");
    assert_eq!(format!("{}", converted_cigar9), "10S80M10S");
}

pub fn read_and_process_overhang(index_dir_adr: &String, bam_file_in: &String,bam_file_out: &String) {
    let mut transcripts_map: HashMap<String, i32> = HashMap::new();
    let mut transcripts: Vec<String> = Vec::new();
    let mut txp_lengths: Vec<i32> = Vec::new();
    let trees = annotation::load_tree(&index_dir_adr, &mut transcripts_map, &mut transcripts, &mut txp_lengths)
            .expect("cannot load the tree!");
    
    let threads = 0;
    let v : Vec<&str> = vec![];
    bam::bam2bam(bam_file_in, bam_file_out, &transcripts, &txp_lengths, &trees, &threads, &200, &15, &v);
}

#[test]
// overhang on one end (either start or end)
fn test_single_exon() {
    
    let index_dir_adr: String = "tests/gencode_v35_index".to_string();
    let sam_file_in = "tests/single_exon.sam".to_string();
    let bam_file_out = "tests/single_exon_toTranscriptome.bam".to_string();
    
    read_and_process_overhang(&index_dir_adr, &sam_file_in, &bam_file_out);

    let mut output_bam = Reader::from_path(&bam_file_out).unwrap();
    let bam_records = output_bam.records();
    let output_header = Reader::from_path(&bam_file_out).unwrap();
    let header_view = output_header.header();

    let pos1 = 0;
    let pos2 = 423;
    let cigar1 = "10S90M";
    let cigar2 = "90M10S";

    for rec in bam_records {
        let record = rec.unwrap();
        let tname_vec = header_view.tid2name(record.tid() as u32).to_vec();
        let tname = std::str::from_utf8(&tname_vec).unwrap();
        match tname {
            "ENST00000612610.4" => {
                assert_eq!(record.pos(), pos1, "Pos1 is wrong! {} {}", record.pos(), pos1);
                assert_eq!(record.cigar().to_string(), cigar1);
            },    
            "ENST00000624081.1" => {
                assert_eq!(record.pos(), pos2, "Pos2 is wrong! {} {}", record.pos(), pos2);
                assert_eq!(record.cigar().to_string(), cigar2);
            },
            _ => continue
        }
    }
}

#[test]
// overhang on both ends
fn test_single_exon2() {
    
    let index_dir_adr: String = "tests/gencode_v35_index".to_string();
    let sam_file_in = "tests/single_exon2.sam".to_string();
    let bam_file_out = "tests/single_exon2_toTranscriptome.bam".to_string();
    
    read_and_process_overhang(&index_dir_adr, &sam_file_in, &bam_file_out);

    let mut output_bam = Reader::from_path(&bam_file_out).unwrap();
    let bam_records = output_bam.records();
    let output_header = Reader::from_path(&bam_file_out).unwrap();
    let header_view = output_header.header();

    let pos1 = 0;
    //let pos2 = 423;
    let cigar1 = "9S76M15S";
    //let cigar2 = "90M10S7

    for rec in bam_records {
        let record = rec.unwrap();
        let tname_vec = header_view.tid2name(record.tid() as u32).to_vec();
        let tname = std::str::from_utf8(&tname_vec).unwrap();
        match tname {  
            "ENST00000624081.1" => {
                assert_eq!(record.pos(), pos1, "Pos1 is wrong! {} {}", record.pos(), pos1);
                assert_eq!(record.cigar().to_string(), cigar1);
            },
            _ => continue
        }
    }
}

#[test]
fn test_multi_exons() {
    
    let index_dir_adr: String = "tests/gencode_v35_index".to_string();
    let sam_file_in = "tests/multi_exons.sam".to_string();
    let bam_file_out = "tests/multi_exons_toTranscriptome.bam".to_string();
    
    read_and_process_overhang(&index_dir_adr, &sam_file_in, &bam_file_out);

    let mut output_bam = Reader::from_path(&bam_file_out).unwrap();
    let bam_records = output_bam.records();
    //let output_header = Reader::from_path(&bam_file_out).unwrap();
    //let header_view = output_header.header();

    let pos1 = 0;
    let pos2 = 0;
    let pos3 = 10;
    let cigar1 = "10S593M10S";
    let cigar2 = "10S583M";
    let cigar3 = "583M10S";

    for rec in bam_records {
        let record = rec.unwrap(); 
        let qname = std::str::from_utf8(&record.qname()).unwrap();

        match qname {  
            "chr21_3exon_10to593to10_ENST00000612610.4" => {
                assert_eq!(record.pos(), pos1, "Pos1 is wrong! {} {}", record.pos(), pos1);
                assert_eq!(record.cigar().to_string(), cigar1);
            },
            "chr21_3exon_10to583_ENST00000612610.4" => {
                assert_eq!(record.pos(), pos2, "Pos2 is wrong! {} {}", record.pos(), pos2);
                assert_eq!(record.cigar().to_string(), cigar2);
            },
            "chr21_3exon_583to10_ENST00000612610.4" => {
                assert_eq!(record.pos(), pos3, "Pos3 is wrong! {} {}", record.pos(), pos3);
                assert_eq!(record.cigar().to_string(), cigar3);
            },
            _ => continue
        }
    }
}

#[test]
// overhang for paired end reads
fn test_paired_reads() {
    
    let index_dir_adr: String = "tests/gencode_v35_index".to_string();
    let sam_file_in = "tests/paired.sam".to_string();
    let bam_file_out = "tests/paired_toTranscriptome.bam".to_string();
    
    read_and_process_overhang(&index_dir_adr, &sam_file_in, &bam_file_out);

    let mut output_bam = Reader::from_path(&bam_file_out).unwrap();
    let bam_records = output_bam.records();
    let output_header = Reader::from_path(&bam_file_out).unwrap();
    let header_view = output_header.header();

    for rec in bam_records {
        let record = rec.unwrap();
        let tname_vec = header_view.tid2name(record.tid() as u32).to_vec();
        let tname = std::str::from_utf8(&tname_vec).unwrap();
        match tname {  
            "ENST00000620481.4" => {
                if record.flags() == 355 {
                    assert_eq!(record.pos(), 0);
                    assert_eq!(record.cigar().to_string(), "10S90M");
                } else {
                    assert_eq!(record.pos(), 110);
                    assert_eq!(record.cigar().to_string(), "91M9S");
                }   
            },
            "ENST00000612610.4" => {
                if record.flags() == 355 {
                    assert_eq!(record.pos(), 0);
                    assert_eq!(record.cigar().to_string(), "10S90M");
                } else {
                    assert_eq!(record.pos(), 110);
                    assert_eq!(record.cigar().to_string(), "91M9S");
                }   
            },
            _ => continue
        }
    }
}