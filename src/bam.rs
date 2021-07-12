use coitrees::{COITree}; //, IntervalNode, SortedQuerent};

use rust_htslib::bam::{Format, Header, Read, Reader, Writer, header}; //, HeaderView};
// use std::collections::HashMap;
// use std::collections::LinkedList;

use crate::annotations;
use annotations::ExonNode;

use crate::intersection;

extern crate fnv;
use fnv::FnvHashMap;



pub fn read_bamfile(input_bam_filename: &String, 
                    // transcripts_map: &HashMap<String, i32>,
                    transcripts: &Vec<String>,
                    txp_lengths: &Vec<i32>,
                    trees: &FnvHashMap::<String, COITree<ExonNode, u32>>) {
    let mut input_bam = Reader::from_path(input_bam_filename).unwrap();
    // let header_ = Header::from_template(input_bam.header());
    let input_bam2 = Reader::from_path(input_bam_filename).unwrap();
    let header_view = input_bam2.header();
    let x = (input_bam).records();
    
    let mut new_header = Header::new();

    // let header_text = String::from_utf8(header_.to_bytes().to_owned()).unwrap().;
    // new_header.push_record(header::HeaderRecord::new(b"HD").push_tag(b"VN", &"1.4"));

    let mut counter = 0;
    for tid in transcripts.iter() {
        let txp_len = txp_lengths[counter];
        new_header.push_record(header::HeaderRecord::new(b"SQ").push_tag(b"SN", &tid).push_tag(b"LN", &txp_len));
        counter += 1;
    }

    let mut output_bam = Writer::from_path("/home/mohsen/test.bam", &new_header, Format::Bam).unwrap();

    // input_bam.set_threads(2).expect("Failed to set number of BAM reading threads to 2.");
    // output_bam.set_threads(5).expect("Failed to set number of BAM writing threads to 4.");

    for r in x {
        let record = r.unwrap();
        let ranges = intersection::find_ranges(&(record.pos() as i32), record.cigar());
        let genome_tname = String::from_utf8(header_view.tid2name(record.tid() as u32).to_vec()).expect("cannot find the tname!");
        if let Some(tree) = trees.get(&genome_tname) {
            let tids = intersection::find_tid(&tree, &ranges);
            if tids.len() > 0 {
                println!("\n {} {} {:?}\n", record.pos(), record.cigar(), genome_tname);
                for tid in tids.iter() {
                    println!("{} {}", tid, transcripts[*tid as usize]);
                }
                println!("done!");
            }
            // record.set_tid(tid);
        }
        output_bam.write(&record).unwrap();
    }
}