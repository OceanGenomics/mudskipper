use coitrees::COITree;

use rust_htslib::bam::{header, record, Format, Header, Read, Reader, Writer};

use crate::annotation;
use annotation::ExonNode;

extern crate bio_types;

use crate::convert;

extern crate fnv;
use fnv::FnvHashMap;

use log::{debug, error, info};

pub fn bam2bam(
    input_bam_filename: &String,
    output_bam_filename: &String,
    transcripts: &Vec<String>,
    txp_lengths: &Vec<i32>,
    trees: &FnvHashMap<String, COITree<ExonNode, u32>>,
    threads_count: &usize,
    max_softlen: &usize,
    required_tags: &Vec<&str>,
) {
    let mut input_bam = Reader::from_path(input_bam_filename).unwrap();
    // let bam_records = input_bam.records();

    // let header_ = Header::from_template(input_bam.header());
    // let header_text = String::from_utf8(header_.to_bytes().to_owned()).unwrap().;

    let input_bam_header = Reader::from_path(input_bam_filename).unwrap();
    let header_view = input_bam_header.header();
    let mut new_header = Header::new();
    // new_header.push_record(header::HeaderRecord::new(b"HD").push_tag(b"VN", &"1.4"));

    info!("Number of reference sequences: {}", transcripts.len());
    let mut counter = 0;
    for tid in transcripts.iter() {
        let txp_len = txp_lengths[counter];
        new_header.push_record(header::HeaderRecord::new(b"SQ").push_tag(b"SN", &tid).push_tag(b"LN", &txp_len));
        counter += 1;
    }

    let mut output_bam = Writer::from_path(output_bam_filename, &new_header, Format::Bam).unwrap();

    if *threads_count >= 2 {
        let threads_count_half = threads_count / 2;
        println!("thread count: {}", threads_count_half);
        input_bam
            .set_threads(threads_count_half)
            .expect("Failed to set number of BAM reading threads.");
        output_bam
            .set_threads(threads_count_half)
            .expect("Failed to set number of BAM writing threads.");
    }

    let mut n = 0;
    let mut mate_wanted = false;
    let mut first_record: record::Record = record::Record::new();
    for rec in input_bam.records() {
        n = n + 1;
        let record = rec.unwrap();
        let qname = String::from_utf8(record.qname().to_vec()).unwrap();
        for tag in required_tags.iter() {
            if record.aux(tag.as_bytes()).is_err() {
                error!("Could not find {} tag for read {}", tag, qname);
                panic!("Some required tags do not exist!");
            }
        }

        debug!("qname: {}", qname);
        if !record.is_paired() || record.is_mate_unmapped() {
            let txp_records = convert::convert_single_end(&record, header_view, transcripts, trees, max_softlen);
            for txp_rec in txp_records.iter() {
                output_bam.write(txp_rec).unwrap();
            }
            mate_wanted = false;
        } else {
            if mate_wanted {
                // this is the second read in pair... PROCESS!
                let txp_records = convert::convert_paired_end(&first_record, &record, header_view, transcripts, txp_lengths, trees, max_softlen);
                for txp_rec in txp_records.iter() {
                    output_bam.write(txp_rec).unwrap();
                }
                mate_wanted = false;
            } else {
                // this is the first read in pair... save and wait for the second read in the pair
                mate_wanted = true;
                first_record = record;
            }
        }
    }
}
