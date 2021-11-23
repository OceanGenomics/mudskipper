use std::convert::TryFrom;

use rust_htslib::bam::{Record, record::Aux, record::CigarString, HeaderView, Read, Reader};

pub struct BAMQueryRecord {
    is_paired: bool,
    first: Vec<Record>,
    second: Vec<Record>,
}

impl BAMQueryRecord {
    pub fn is_paired(&self) -> bool {
        self.is_paired
    }

    pub fn get_first(&self) -> &Vec<Record> {
        &self.first
    }

    pub fn get_second(&self) -> &Vec<Record> {
        &self.second
    }
}

pub struct BAMQueryRecordReader {
    bam_reader: Reader,
    header: HeaderView,
    last_qname: String,
    record_list: Vec<Record>,
    supp_list: Vec<Record>,
}

impl BAMQueryRecordReader {
    pub fn new(bam_filename: &String) -> BAMQueryRecordReader {
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

    // Returns the list of all chimeric alignments in the SA tag (by convention the first one is the primary alignment)
    #[allow(dead_code)]
    fn get_records_from_sa_tag(&self, sa_tag: &str) -> Vec<Record> {
        let mut sa_records: Vec<Record> = Vec::new();
        for aln_str in sa_tag.split(";") {
            if aln_str.is_empty() == false {
                let mut brecord = Record::new();
                let tag_vec: Vec<&str> = aln_str.split(",").collect();

                let cigar: CigarString = CigarString::try_from(tag_vec[3].as_bytes()).expect("Unable to parse cigar string.");
                brecord.set("".as_bytes(), Some(&cigar), "".as_bytes(), "".as_bytes());
                brecord.set_tid(self.header.tid(tag_vec[0].as_bytes()).expect("Cannot find tid for SA alignment!") as i32);
                brecord.set_pos(tag_vec[1].parse::<i64>().expect("Cannot parse position for SA alignment!") -1);
                if tag_vec[2].eq("-") == true {
                    brecord.set_reverse();
                } else {
                    brecord.unset_reverse();
                }
                brecord.set_mapq(tag_vec[4].parse::<u8>().expect("Cannot parse MAPQ for SA alignment!"));
                sa_records.push(brecord);
            }
        }
        sa_records
    }

    // By convention, the first alignment in the SA tag refers to the primary alignment
    fn get_primary_record_of_sa_tag(&self, sa_tag: &str) -> Record {
        let mut brecord = Record::new();
        for aln_str in sa_tag.split(";") {
            if aln_str.is_empty() == false {
                let tag_vec: Vec<&str> = aln_str.split(",").collect();

                let cigar: CigarString = CigarString::try_from(tag_vec[3].as_bytes()).expect("Unable to parse cigar string.");
                brecord.set("".as_bytes(), Some(&cigar), "".as_bytes(), "".as_bytes());
                brecord.set_tid(self.header.tid(tag_vec[0].as_bytes()).expect("Cannot find tid for SA alignment!") as i32);
                brecord.set_pos(tag_vec[1].parse::<i64>().expect("Cannot parse position for SA alignment!") -1);
                if tag_vec[2].eq("-") == true {
                    brecord.set_reverse();
                } else {
                    brecord.unset_reverse();
                }
                brecord.set_mapq(tag_vec[4].parse::<u8>().expect("Cannot parse MAPQ for SA alignment!"));
                return brecord
            }
        }
        brecord
    }

    fn records_equal(a: &Record, b: &Record) -> bool {
        // println!("\tCOMPARE {}:{}:{}", String::from_utf8(self.header.tid2name(a.tid() as u32).to_vec()).expect("cannot find the tname!"), a.pos(), a.is_reverse());
        println!("\tCOMPARE {}:{} {}:{} cigar:{}", a.pos(), a.is_reverse(), b.pos(), b.is_reverse(), a.cigar() == b.cigar());
        a.tid() == b.tid()
        && a.pos() == b.pos()
        && a.is_reverse() == b.is_reverse()
        && a.cigar() == b.cigar()
    }

    fn group_records(&self) -> Vec<BAMQueryRecord> {
        let mut record_groups: Vec<BAMQueryRecord> = Vec::new();
        // find primary alignment of each supplementary alignment
        let mut primary_of_supp: Vec<Record> = Vec::new();
        for supp in self.supp_list.iter() {
            if let Ok(Aux::String(sa_tag)) = supp.aux(b"SA") {
                primary_of_supp.push(self.get_primary_record_of_sa_tag(sa_tag))
            }
            else {
                panic!("Error reading SA tag for query {}", String::from_utf8(supp.qname().to_vec()).expect("cannot find the qname!"));
            }
        }
        // 
        let mut supp_assigned: Vec<bool> = vec![false; self.supp_list.len()];
        let mut i = 0;
        while i < self.record_list.len() {
            let mut first_vec: Vec<Record> = vec![];
            let mut second_vec: Vec<Record> = vec![];
            if self.record_list[i].is_paired() == false { // single-end
                // add the primary alignment first
                first_vec.push(self.record_list[i].to_owned());
                // println!("INSIDE\tSE {}:{}", 
                //     String::from_utf8(self.header.tid2name(self.record_list[i].tid() as u32).to_vec()).expect("cannot find the tname!"), 
                //     self.record_list[i].pos()
                // );
                // search for possible matching of a supplementary alignment with the alignment in first_vec
                for j in 0..supp_assigned.len() {
                    if supp_assigned[j] == false && BAMQueryRecordReader::records_equal(&self.record_list[i], &primary_of_supp[j]) {
                        println!("ADDING to first_vec");
                        first_vec.push(self.supp_list[j].to_owned());
                        supp_assigned[j] = true;
                    }
                }
                // for supp in self.supp_list.iter() {
                //     if let Ok(Aux::String(sa_tag)) = supp.aux(b"SA") {
                //         println!("SUPP\t{}", sa_tag);
                        
                //     }
                //     else {
                //         panic!("Error reading SA tag for query {}", String::from_utf8(supp.qname().to_vec()).expect("cannot find the qname!"));
                //     }
                // }
                record_groups.push(BAMQueryRecord {
                    is_paired: false,
                    first: first_vec,
                    second: second_vec,
                });
                i += 1;
            } else { // paired-end
                // add the primary alignments first
                first_vec.push(self.record_list[i].to_owned());
                second_vec.push(self.record_list[i+1].to_owned());
                // println!("INSIDE\tPE {}:{} {}:{}", 
                //     String::from_utf8(self.header.tid2name(self.record_list[i].tid() as u32).to_vec()).expect("cannot find the tname!"),
                //     self.record_list[i].pos(), 
                //     String::from_utf8(self.header.tid2name(self.record_list[i+1].tid() as u32).to_vec()).expect("cannot find the tname!"),
                //     self.record_list[i+1].pos(), 
                // );
                // search for possible matching of a supplementary alignment with the alignment in first_vec or second_vec
                for j in 0..supp_assigned.len() {
                    if supp_assigned[j] == false {
                        if BAMQueryRecordReader::records_equal(&self.record_list[i], &primary_of_supp[j]) {
                            println!("ADDING to first_vec");
                            first_vec.push(self.supp_list[j].to_owned());
                            supp_assigned[j] = true;
                        } else if BAMQueryRecordReader::records_equal(&self.record_list[i+1], &primary_of_supp[j]) {
                            println!("ADDING to second_vec");
                            second_vec.push(self.supp_list[j].to_owned());
                            supp_assigned[j] = true;
                        }
                    }
                }
                // // search for possible matching supplementary alignment
                // for supp in self.supp_list.iter() {
                //     match supp.aux(b"SA") {
                //         Ok(value) => {
                //             if let Aux::String(v) = value {
                //                 println!("SUPP\t{}", v);
                //             }
                //         }
                //         Err(e) => {
                //             panic!("Error reading SA tag: {}", e);
                //         }
                //     }
                // }
                record_groups.push(BAMQueryRecord {
                    is_paired: true,
                    first: first_vec,
                    second: second_vec,
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