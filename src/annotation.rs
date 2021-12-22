use bio::io::gff;
use coitrees::{COITree, IntervalNode};
use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::time::Instant;

extern crate bio_types;
extern crate fnv;

use bio_types::strand::Strand;
use fnv::FnvHashMap;

use indicatif::ProgressBar;
use linecount::count_lines;
use log;

pub fn read(ann_file_adr: &String) -> Result<gff::Reader<File>, Box<dyn Error>> {
    let ann_file_adr_split: Vec<&str> = ann_file_adr.split(".").collect();
    let file_type: &str = ann_file_adr_split.last().copied().unwrap_or("default string");
    log::info!("reading the {} file and building the tree.", file_type);
    let ann_type: gff::GffType = if file_type == "gtf" {
        gff::GffType::GTF2
    } else {
        if file_type == "gff3" || file_type == "gff" {
            gff::GffType::GFF3
        } else {
            gff::GffType::GFF2
        }
    };

    return Ok(gff::Reader::from_file(ann_file_adr, ann_type).expect("Error in reading annotation file."));
}

pub struct ExonNode {
    pub start: i32,
    pub end: i32,
    pub tid: i32,
    pub tpos_start: i32,
    pub strand: Strand,
}

impl std::fmt::Display for ExonNode {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "(start: {}, end : {}, tid : {}, tpos: {}, strand: {})",
            self.start, self.end, self.tid, self.tpos_start, self.strand
        )
    }
}

impl Clone for ExonNode {
    fn clone(&self) -> Self {
        let new_exon: ExonNode = ExonNode {
            start: self.start,
            end: self.end,
            tid: self.tid,
            tpos_start: self.tpos_start,
            strand: self.strand,
        };
        return new_exon;
    }
}

pub fn load_tree(
    transcripts_map: &mut HashMap<String, i32>,
    transcripts: &mut Vec<String>,
    txp_lengths: &mut Vec<i32>,
) -> Result<FnvHashMap<String, COITree<ExonNode, u32>>, Box<dyn Error>> {
    log::info!("Loading parsed GTF...");
    // load info
    {
        let mut ifile = File::open("parsed_gtf.map").unwrap();
        let mut ireader = BufReader::new(ifile);
        for line in ireader.lines() {
            let line_str = line.unwrap();
            let cols: Vec<&str> = line_str.trim().split("\t").collect();
            transcripts_map.insert(cols[0].to_string(), cols[1].parse::<i32>().unwrap());
        }
        ifile = File::open("parsed_gtf.name").unwrap();
        ireader = BufReader::new(ifile);
        for line in ireader.lines() {
            transcripts.push(line.unwrap().trim().to_string());
        }
        ifile = File::open("parsed_gtf.len").unwrap();
        ireader = BufReader::new(ifile);
        for line in ireader.lines() {
            txp_lengths.push(line.unwrap().parse::<i32>().unwrap());
        }
    }
    //
    let mut trees = FnvHashMap::<String, COITree<ExonNode, u32>>::default();
    let ifile = File::open("parsed_gtf.exon").unwrap();
    let ireader = BufReader::new(ifile);
    let mut last_name: String = String::from("");
    let mut node_vec: Vec<IntervalNode<ExonNode, u32>> = Vec::new();
    for line in ireader.lines() {
        let line_str = line.unwrap();
        let cols: Vec<&str> = line_str.trim().split("\t").collect();
        let seq_name = cols[0];
        let exon_start: i32 = cols[1].parse().unwrap();
        let exon_end: i32 = cols[2].parse().unwrap();
        let exon_tid: i32 = cols[3].parse().unwrap();
        let exon_tpos: i32 = cols[4].parse().unwrap();
        let exon_strand: Strand = cols[5].parse().unwrap();
        let exon: ExonNode = ExonNode {
            start: exon_start,
            end: exon_end,
            tid: exon_tid,
            tpos_start: exon_tpos,
            strand: exon_strand,
        };

        if seq_name != last_name {
            // add exisiting vector to the tree and initiate a new one
            if node_vec.len() > 0 {
                trees.insert(last_name, COITree::new(node_vec.clone()));
            }
            //
            last_name = seq_name.to_string();
            node_vec.clear();
            node_vec.push(IntervalNode::new(exon_start, exon_end, exon))
        } else {
            // add to the vector
            node_vec.push(IntervalNode::new(exon_start, exon_end, exon))
        }
    }
    // add last vector to the tree
    if node_vec.len() > 0 {
        trees.insert(last_name, COITree::new(node_vec.clone()));
    }
    return Ok(trees);
}

pub fn build_tree(
    ann_file_adr: &String,
    transcripts_map: &mut HashMap<String, i32>,
    transcripts: &mut Vec<String>,
    txp_lengths: &mut Vec<i32>,
) -> Result<FnvHashMap<String, COITree<ExonNode, u32>>, Box<dyn Error>> {
    let mut nodes = FnvHashMap::<String, Vec<IntervalNode<ExonNode, u32>>>::default();
    let a = Instant::now();
    let reader = read(ann_file_adr);
    let mut tid: i32 = 0;
    let mut tpos: i32 = 0;

    // let gtf_records_count = read(ann_file_adr).expect("Error reading file.")
    //                                          .records().count();
    let gtf_records_count = count_lines(File::open(ann_file_adr).unwrap()).unwrap();
    let pb = ProgressBar::new(gtf_records_count as u64);
    for record in reader.expect("Error reading file.").records() {
        pb.inc(1);
        let rec = record.ok().expect("Error reading record.");
        let features = rec.attributes();
        let tname_key: String = "transcript_id".to_string();
        if rec.feature_type() == "exon" && features.contains_key(&tname_key) {
            if (features.contains_key("exon_number") && features["exon_number"] == "1") || (features.contains_key("exon") && features["exon"] == "1")
            {
                tpos = 0;
            }
            let seqname = rec.seqname().to_string();
            // log::debug!("{:?}\t{:?}", rec.feature_type(), seqname);
            let exon_start = (*rec.start() - 1) as i32;
            let exon_end = (*rec.end() - 1) as i32;
            let exon_len = exon_end - exon_start + 1;
            let exon_strand = rec.strand(); //.unwrap();
            let exon_strand = match exon_strand {
                Some(strand) => strand,
                None => {
                    log::debug!("The gtf/gff record doesn't specify the strand, will be ignored.");
                    continue;
                }
            };
            let features = rec.attributes();
            let tname = &features[&tname_key];
            if features.contains_key(&tname_key) {
                if !transcripts_map.contains_key(&tname.to_string()) {
                    transcripts_map.insert(tname.to_string(), tid);
                    transcripts.push(tname.to_string());
                    txp_lengths.push(exon_len);
                    tid += 1;
                } else {
                    let _tid = transcripts_map[&tname.to_string()] as usize;
                    txp_lengths[_tid] += exon_len;
                }
            } else {
            }
            let exon: ExonNode = ExonNode {
                start: exon_start,
                end: exon_end,
                tid: transcripts_map[&tname.to_string()],
                tpos_start: tpos,
                strand: exon_strand,
            };
            let node_arr = if let Some(node_arr) = nodes.get_mut(&seqname[..]) {
                node_arr
            } else {
                nodes.entry(seqname).or_insert(Vec::new())
            };
            node_arr.push(IntervalNode::new(exon_start, exon_end, exon));

            // Update the tpos for the next exon
            // println!("start:{} end:{} exon_number:{}  tpos:{}",
            // exon_start, exon_end, features["exon_number"], tpos);
            tpos += exon_len;
        }
    }
    pb.finish_with_message("finish reading the file");

    // save info
    {
        let mut outfile = File::create("parsed_gtf.map").unwrap();
        for (seqname, seqid) in transcripts_map.iter() {
            write!(outfile, "{}\t{}\n", seqname, seqid).unwrap();
        }
        outfile = File::create("parsed_gtf.name").unwrap();
        for tname in transcripts.iter() {
            write!(outfile, "{}\n", tname).unwrap();
        }
        outfile = File::create("parsed_gtf.len").unwrap();
        for tlen in txp_lengths.iter() {
            write!(outfile, "{}\n", tlen).unwrap();
        }
    }

    log::info!("building the tree");
    let mut outfile = File::create("parsed_gtf.exon").unwrap();
    let mut trees = FnvHashMap::<String, COITree<ExonNode, u32>>::default();
    for (seqname, seqname_nodes) in nodes {
        // write!(outfile, "{}\t{}\n", seqname, seqname_nodes.len()).unwrap();
        for inode in seqname_nodes.iter() {
            write!(
                outfile,
                "{}\t{}\t{}\t{}\t{}\t{}\n",
                seqname, inode.metadata.start, inode.metadata.end, inode.metadata.tid, inode.metadata.tpos_start, inode.metadata.strand
            )
            .unwrap();
        }
        trees.insert(seqname, COITree::new(seqname_nodes));
    }
    let b = Instant::now();
    log::info!("Time to build the tree: {:?}", b - a);
    return Ok(trees);
}
