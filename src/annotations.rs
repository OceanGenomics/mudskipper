use std::time::Instant;
use std::error::Error;
use bio::io::gff;
use coitrees::{COITree, IntervalNode};
use std::collections::HashMap;

extern crate bio_types;
use bio_types::strand::Strand;

extern crate fnv;
use fnv::FnvHashMap;

use log::{info, debug};
use indicatif::ProgressBar;
use linecount::count_lines;

type GenericError = Box<dyn Error>;

pub fn read(ann_file_adr: &String) -> Result<gff::Reader<std::fs::File>, GenericError> {
    let ann_file_adr_split: Vec<&str> = ann_file_adr.split(".").collect();
    let file_type: &str = ann_file_adr_split.last().copied().unwrap_or("default string");
    info!("reading the {} file and building the tree.", file_type);
    let ann_type: gff::GffType = if file_type == "gtf" { gff::GffType::GTF2 } 
                                 else { if file_type == "gff3" || file_type == "gff" { gff::GffType::GFF3 }
                                        else { gff::GffType::GFF2 } };

    return Ok(gff::Reader::from_file(ann_file_adr, ann_type).expect("Error in reading annotation file."));
}

pub struct ExonNode {
    pub start: i32,
    pub end: i32,
    pub tid: i32,
    pub tpos_start: i32,
    pub strand: Strand
}

impl std::fmt::Display for ExonNode {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "(start: {}, end : {}, tid : {}, tpos: {}, strand: {})", 
                    self.start, self.end, self.tid, self.tpos_start, self.strand)
    }
}

impl Clone for ExonNode {
    fn clone(&self) -> Self {
        let new_exon: ExonNode = ExonNode{start: self.start, 
                                          end: self.end,
                                          tid: self.tid,
                                          tpos_start: self.tpos_start,
                                          strand: self.strand 
        };
        return new_exon;
    }
}

pub fn build_tree(ann_file_adr: &String, 
                transcripts_map: &mut HashMap<String, i32>,
                transcripts: &mut Vec<String>,
                txp_lengths: &mut Vec<i32>) 
    -> Result<FnvHashMap<String, COITree<ExonNode, u32>>, GenericError> {
    
    let mut nodes = FnvHashMap::<String, Vec<IntervalNode<ExonNode, u32>>>::default();
    let a = Instant::now();
    let reader = read(ann_file_adr);
    let mut tid: i32 = 0;
    let mut tpos: i32 = 0;

    // let gtf_records_count = read(ann_file_adr).expect("Error reading file.")
    //                                          .records().count();
    let gtf_records_count = count_lines(std::fs::File::open(ann_file_adr).unwrap()).unwrap();
    let pb = ProgressBar::new(gtf_records_count as u64);
    for record in reader.expect("Error reading file.").records() {
        pb.inc(1);
        let rec = record.ok().expect("Error reading record.");
        let features = rec.attributes();
        let tname_key : String = "transcript_id".to_string();
        if rec.feature_type() == "exon" && features.contains_key(&tname_key) {
            if (features.contains_key("exon_number") && features["exon_number"] == "1") ||
                (features.contains_key("exon") && features["exon"] == "1") {
                tpos = 0;
            }
            let seqname = rec.seqname().to_string();
            debug!("\r{:?}\t{:?}", rec.feature_type(), seqname);
            let exon_start = *rec.start() as i32;
            let exon_end = *rec.end() as i32;
            let exon_len = exon_end-exon_start+1;
            let exon_strand = rec.strand();//.unwrap();
            let exon_strand = match exon_strand {
                Some(strand) => strand,
                None => {
                    debug!("The gtf/gff record doesn't specify the strand, will be ignored.");
                    continue;
                },
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
            let exon: ExonNode = ExonNode{start: exon_start,
                                            end: exon_end,
                                            tid: transcripts_map[&tname.to_string()],
                                            tpos_start: tpos,
                                            strand: exon_strand};
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

    info!("building the tree");
    let mut trees = FnvHashMap::<String, COITree<ExonNode, u32>>::default();    
    for (seqname, seqname_nodes) in nodes {
        trees.insert(seqname, COITree::new(seqname_nodes));
    }
    let b = Instant::now();
    info!("Time to build the tree: {:?}", b-a);
    return Ok(trees);
}