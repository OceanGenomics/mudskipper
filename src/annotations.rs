use std::time::Instant;
use std::error::Error;
use bio::io::gff;
use coitrees::{COITree, IntervalNode}; //, SortedQuerent};
use std::collections::HashMap;
// use std::collections::LinkedList;

extern crate fnv;
use fnv::FnvHashMap;

type GenericError = Box<dyn Error>;

pub fn read(ann_file_adr: &String) -> Result<gff::Reader<std::fs::File>, GenericError > {
    let ann_file_adr_split: Vec<&str> = ann_file_adr.split(".").collect();
    let file_type: &str = ann_file_adr_split.last().copied().unwrap_or("default string");
    println!("{}", file_type);
    let ann_type: gff::GffType = if file_type == "gtf" { gff::GffType::GTF2 } 
                                 else { if file_type == "gff3" { gff::GffType::GFF3 }
                                        else { gff::GffType::GFF2 } };

    return Ok(gff::Reader::from_file(ann_file_adr, ann_type).expect("Error reading file."))
}

pub struct ExonNode {
    // transcript: String,
    start: i32,
    end: i32,
    pub tid: i32
}

impl std::fmt::Display for ExonNode {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "(start: {}, end : {}, tid : {})", self.start, self.end, self.tid)
    }
}

/*impl Copy for ExonNode {

}*/

impl Clone for ExonNode {
    fn clone(&self) -> Self {
        let new_exon: ExonNode = ExonNode{start: self.start, 
                                            end: self.end,
                                            tid: self.tid};
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
    let mut n: i32 = 0;
    let mut tid: i32 = 0;
    for record in reader.expect("Error reading file.").records() {
        let rec = record.ok().expect("Error reading record.");
        if rec.feature_type() == "exon" {
            n += 1;
            let seqname = rec.seqname().to_string();
            print!("\r{}\t{:?}\t{:?}", n, rec.feature_type(), seqname);
            let exon_start = *rec.start() as i32;
            let exon_end = *rec.end() as i32;
            let exon_len = exon_end-exon_start+1;
            let features = rec.attributes();
            let tname_key : String = "transcript_id".to_string();
            let tname = &features[&tname_key];
            if features.contains_key(&tname_key) {                
                // println!("{}", tname);
               if !transcripts_map.contains_key(&tname.to_string()) {
                    transcripts_map.insert(tname.to_string(), tid);
                    transcripts.push(tname.to_string());                    
                    txp_lengths.push(exon_len);
                    tid += 1;
                } else {
                    let _tid = transcripts_map[&tname.to_string()] as usize;
                    txp_lengths[_tid] += exon_len;
                }
            }   
            let exon: ExonNode = ExonNode{start: exon_start,
                                            end: exon_end,
                                            tid: transcripts_map[&tname.to_string()]};
            let node_arr = if let Some(node_arr) = nodes.get_mut(&seqname[..]) {
                node_arr
            } else {
                nodes.entry(seqname).or_insert(Vec::new())
            };
            node_arr.push(IntervalNode::new(exon_start, exon_end, exon));
        }
        if n > 1000 {
            break;
        }
    }
    println!("\nn: {}", n);
    let mut trees = FnvHashMap::<String, COITree<ExonNode, u32>>::default();
    for (seqname, seqname_nodes) in nodes {
        trees.insert(seqname, COITree::new(seqname_nodes));
    }
    let b = Instant::now();
    println!("Time to build: {:?}", b-a);
    return Ok(trees);
}

pub fn test_tree(ann_file_adr: &String, trees: &FnvHashMap::<String, COITree<ExonNode, u32>>) {
    let reader = read(ann_file_adr);
    let mut n: i32 = 0;
    let a = Instant::now();
    for record in reader.expect("Error reading file.").records() {
        let rec = record.ok().expect("Error reading record.");
        let seqname = rec.seqname().to_string();
        if rec.feature_type() == "exon" {
            n += 1;
            let exon_start = *rec.start() as i32;
            let exon_end = *rec.end() as i32;
            if let Some(seqname_tree) = trees.get(&seqname) {
                let countcov = seqname_tree.coverage(exon_start, exon_end);
                let count = countcov.0;
                let cov = countcov.1;
                if count == 0 || cov == 0 {
                    println!("Error happend!")
                } else {
                    // println!("{} {}", count, cov);
                }
            }
        }
        if n > 100 {
        //    break;
        }
    }
    let b = Instant::now();
    println!("Time to query: {:?}", b-a);
}