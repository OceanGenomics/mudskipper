use std::io::{self, Write};

use coitrees::COITree;

use rust_htslib::bam::{header, Record, Format, Header, HeaderView, Read, Reader, Writer};

use crate::annotation;
use annotation::ExonNode;

extern crate bio_types;

use crate::convert;

extern crate fnv;
use fnv::FnvHashMap;

use log::{debug, error, info};

struct BAMQueryRecord {
    is_paired: bool,
    first: Vec<Record>,
    second: Vec<Record>,
}

struct BAMQueryRecordReader {
    bam_reader: Reader,
    header: HeaderView,
    last_qname: String,
    record_list: Vec<Record>,
    supp_list: Vec<Record>,
}

impl BAMQueryRecordReader {
    fn new(bam_filename: &String) -> BAMQueryRecordReader {
        let mut breader = Reader::from_path(bam_filename).unwrap();
        let hv = breader.header().to_owned();
        // init record_list by finding the first non-supplementary record
        let mut last_qname = String::from("");
        let mut r_list: Vec<Record> = Vec::new();
        let mut s_list: Vec<Record> = Vec::new();
        let mut brecord = Record::new();
        while let Some(res) = breader.read(&mut brecord) {
            res.expect("Failed to parse BAM record");
            last_qname = String::from_utf8(brecord.qname().to_vec()).unwrap();
            if brecord.flags() & 0x800 != 0 {
                s_list.push(brecord.to_owned());
            } else {
                r_list.push(brecord.to_owned());
                break; // break as soon as the first record has been seen
            }
        }
        if r_list.is_empty() {
            panic!("Could not find any non-supplementary records in the BAM file!");
        }
        BAMQueryRecordReader {
            bam_reader: breader,
            header: hv,
            last_qname: last_qname,
            record_list: r_list,
            supp_list: s_list,
        }
    }

    pub fn get_header(&self) -> &HeaderView {
        &self.header
    }

    fn group_records(&self) -> Vec<BAMQueryRecord> {
        let mut record_groups: Vec<BAMQueryRecord> = Vec::new();
        // let mut mate_wanted = false;
        // let mut first_record: Record;
        // let mut second_record: Record;
        let mut i = 0;
        while i < self.record_list.len() {
            if self.record_list[i].is_paired() == false {
                // single-end
                record_groups.push(BAMQueryRecord {
                    is_paired: false,
                    first: vec![self.record_list[i].to_owned()],
                    second: vec![],
                });
                // println!("INSIDE\tSE {}:{}", 
                //     String::from_utf8(self.header.tid2name(self.record_list[i].tid() as u32).to_vec()).expect("cannot find the tname!"), 
                //     self.record_list[i].pos()
                // );
                i += 1;
            } else {
                // paired-end
                record_groups.push(BAMQueryRecord {
                    is_paired: true,
                    first: vec![self.record_list[i].to_owned()],
                    second: vec![self.record_list[i+1].to_owned()],
                });
                // println!("INSIDE\tPE {}:{} {}:{}", 
                //     String::from_utf8(self.header.tid2name(self.record_list[i].tid() as u32).to_vec()).expect("cannot find the tname!"),
                //     self.record_list[i].pos(), 
                //     String::from_utf8(self.header.tid2name(self.record_list[i+1].tid() as u32).to_vec()).expect("cannot find the tname!"),
                //     self.record_list[i+1].pos(), 
                // );
                i += 2;
            }
        }
        record_groups
    }

    pub fn get_next_query_records(&mut self) -> Vec<BAMQueryRecord> {
        let mut brecord = Record::new();
        let query_records: Vec<BAMQueryRecord>;
        while let Some(res) = self.bam_reader.read(&mut brecord) {
            res.expect("Failed to parse BAM record");
            let qname = String::from_utf8(brecord.qname().to_vec()).unwrap();
            if qname != self.last_qname { // a new query
                // process the last group
                // println!("INSIDE\tqname: {}", self.last_qname);
                // for r in self.record_list.iter() {
                //     println!("INSIDE\t\t{}:{}", String::from_utf8(self.header.tid2name(r.tid() as u32).to_vec()).expect("cannot find the tname!"), r.pos());
                // }
                // for r in self.supp_list.iter() {
                //     println!("INSIDE\t\tsupp\t{}:{}", String::from_utf8(self.header.tid2name(r.tid() as u32).to_vec()).expect("cannot find the tname!"), r.pos());
                // }
                query_records = self.group_records();
                // reset
                self.last_qname = qname;
                self.record_list.clear();
                self.supp_list.clear();
                // add the current record
                if brecord.flags() & 0x800 != 0 {
                    self.supp_list.push(brecord.to_owned());
                } else {
                    self.record_list.push(brecord.to_owned());
                }
                return query_records;
            } else { // same query
                if brecord.flags() & 0x800 != 0 {
                    self.supp_list.push(brecord.to_owned());
                } else {
                    self.record_list.push(brecord.to_owned());
                }
            }
        }
        // process the last record
        // println!("INSIDE\tqname: {}", self.last_qname);
        // for r in self.record_list.iter() {
        //     println!("INSIDE\t\t{}:{}", String::from_utf8(self.header.tid2name(r.tid() as u32).to_vec()).expect("cannot find the tname!"), r.pos());
        // }
        // for r in self.supp_list.iter() {
        //     println!("INSIDE\t\tsupp\t{}:{}", String::from_utf8(self.header.tid2name(r.tid() as u32).to_vec()).expect("cannot find the tname!"), r.pos());
        // }
        query_records = self.group_records();
        // reset
        self.record_list.clear();
        self.supp_list.clear();
        query_records
    }
}

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
    let mut br = BAMQueryRecordReader::new(input_bam_filename);
    let hv = br.get_header().to_owned();
    loop {
        let ret_vec = br.get_next_query_records();
        if ret_vec.is_empty() {
            break;
        }
        // println!("qname: {} {}", r.first[0].qname(), r.second[0].qname());
        for r in ret_vec.iter() {
            if r.is_paired {
                println!("\tPE  query {} {}", String::from_utf8(r.first[0].qname().to_vec()).expect("cannot find the qname!"), String::from_utf8(r.second[0].qname().to_vec()).expect("cannot find the qname!"));
                // println!("\tPE target {}:{} {}:{}", String::from_utf8(hv.tid2name(r.first[0].tid() as u32).to_vec()).expect("cannot find the tname!"), r.first[0].pos(), String::from_utf8(hv.tid2name(r.second[0].tid() as u32).to_vec()).expect("cannot find the tname!"), r.second[0].pos());
                print!("\tPE target");
                if r.first[0].is_unmapped() == false {
                    print!(" {}:{}", String::from_utf8(hv.tid2name(r.first[0].tid() as u32).to_vec()).expect("cannot find the tname!"), r.first[0].pos());
                } else {
                    print!(" {}:{}", r.first[0].tid(), r.first[0].pos());
                }
                if r.second[0].is_unmapped() == false {
                    print!(" {}:{}", String::from_utf8(hv.tid2name(r.second[0].tid() as u32).to_vec()).expect("cannot find the tname!"), r.second[0].pos());
                }
                print!("\n");
                io::stdout().flush().unwrap();
            } else {
                println!("\tSE  query {}", String::from_utf8(r.first[0].qname().to_vec()).expect("alaki"));
                if r.first[0].is_unmapped() {
                    println!("\tSE target {}:{}", String::from_utf8(hv.tid2name(r.first[0].tid() as u32).to_vec()).expect("cannot find the tname!"), r.first[0].pos());
                }
            }
        }
        println!("")
    }
}

pub fn bam2bam_2(
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

    // let input_bam_header = Reader::from_path(input_bam_filename).unwrap();
    let header_view = input_bam.header().to_owned();
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
    let mut first_record: Record = Record::new();
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
            let txp_records = convert::convert_single_end(&record, &header_view, transcripts, trees, max_softlen);
            for txp_rec in txp_records.iter() {
                output_bam.write(txp_rec).unwrap();
            }
            mate_wanted = false;
        } else {
            if mate_wanted {
                // this is the second read in pair... PROCESS!
                let txp_records = convert::convert_paired_end(&first_record, &record, &header_view, transcripts, txp_lengths, trees, max_softlen);
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
