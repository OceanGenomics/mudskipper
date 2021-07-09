use rust_htslib::bam::{Format, Header, Read, Reader, Writer, header};
use std::collections::HashMap;
use std::collections::LinkedList;

pub fn find_tid() -> i32 {
    return 0;
}

pub fn read_bamfile(input_bam_filename: &String, 
                    transcripts_map: &HashMap<String, i32>,
                    transcripts: &LinkedList<String>,
                    txp_lengths: &Vec<i32>) {
    let mut input_bam = Reader::from_path(input_bam_filename).unwrap();
    let header_ = Header::from_template(input_bam.header());
    
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

    for r in input_bam.records() {
        let mut record = r.unwrap();

        let tid = find_tid();
        record.set_tid(tid);
        // println!("{}", record.name());

        output_bam.write(&record).unwrap();
    }
}