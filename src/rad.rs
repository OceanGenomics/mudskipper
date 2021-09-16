// use std::io::{stdout, BufReader, BufWriter, Cursor, Seek, SeekFrom, Write};

use crate::annotations;

extern crate fnv;

use std::fs;
use std::fs::File;
// use std::path::Path;
use log::{info, error};

use annotations::ExonNode;

use std::io::{BufWriter, Cursor, Write, Seek, SeekFrom};
use rust_htslib::{bam, bam::record::Aux, bam::Read};
use fnv::FnvHashMap;
use coitrees::{COITree};

// struct RadRecord {
//     pub fw: bool,
//     pub qname: String,
//     pub tid_list: Vec<i32>
// }

fn bam_peek(bam_path: &String) -> bam::Record {
    let mut bam = bam::Reader::from_path(&bam_path).unwrap();
    bam.set_threads(1).unwrap();
    // 
    let mut rec = bam::Record::new();
    let first_record_exists = bam.read(&mut rec).is_some();
    if !first_record_exists {
        error!("bam file had no records!");
        std::process::exit(1);
    }
    rec
}

pub fn bam2rad(input_bam_filename: &String, 
                 output_rad_filename: &String,
                 transcripts: &Vec<String>,
                 txp_lengths: &Vec<i32>,
                 trees: &FnvHashMap::<String, COITree<ExonNode, u32>>,
                 threads_count: &usize,
                 max_softlen: &usize) {

    let ofile = File::create(&output_rad_filename).unwrap();

    // {
    //     let bam_bytes = fs::metadata(&input_bam_filename).unwrap().len();
    //     info!("Bam file size in bytes {:?}", bam_bytes);
    // }

    // file writer
    let mut owriter = BufWriter::with_capacity(1048576, ofile);
    // intermediate buffer
    let mut data = Cursor::new(vec![]);
    let first_rec = bam_peek(input_bam_filename);

    ///////////////////////////////////////// header section
    {
        // is paired-end
        let is_paired: u8;
        if first_rec.is_paired() {
            is_paired = 1;
        } else {
            is_paired = 0;
        }
        data.write_all(&is_paired.to_le_bytes())
            .expect("couldn't write to output file");
        // number of references
        let ref_count = transcripts.len() as u64;
        println!("ref_count: {}", ref_count);
        data.write_all(&ref_count.to_le_bytes())
            .expect("couldn't write to output file");
        // list of reference names
        for tname in transcripts.iter() {            
            let name_len = tname.len() as u16;
            data.write_all(&name_len.to_le_bytes())
                .expect("coudn't write to output file");
            data.write_all(tname.as_bytes())
                .expect("coudn't write to output file");
        }
        // number of chunks
        let initial_num_chunks = 0u64;
        data.write_all(&initial_num_chunks.to_le_bytes())
            .expect("coudn't write to output file");
    }

    // keep a pointer to header pos
    let end_header_pos =
    data.seek(SeekFrom::Current(0)).unwrap() - std::mem::size_of::<u64>() as u64;
    // check header position
    info!("end header pos: {:?}", end_header_pos);

    ///////////////////////////////////////// tag definition section
    {
        // FILE-LEVEL tags
        // TODO: assumes bulk sample only for now
        let mut num_tags = 0u16; // no file-level tags for bulk
        data.write_all(&num_tags.to_le_bytes())
            .expect("coudn't write to output file");

        // READ-LEVEL tags
        // TODO: assumes bulk sample only for now
        num_tags = 0u16; // no read-level tags for bulk
        data.write_all(&num_tags.to_le_bytes())
            .expect("coudn't write to output file");

        // ALIGNMENT-LEVEL tags
        // TODO: assumes bulk sample only for now
        num_tags = 0u16;
        data.write_all(&num_tags.to_le_bytes())
            .expect("couldn't write to output file");
    }

    ///////////////////////////////////////// file-level tag values
    {
        // nothing for bulk!
    }

    let alaki = 255u8;
    data.write_all(&alaki.to_le_bytes())
    .expect("coudn't write to output file");
    data.write_all(&alaki.to_le_bytes())
    .expect("coudn't write to output file");
    data.write_all(&alaki.to_le_bytes())
    .expect("coudn't write to output file");
    data.write_all(&alaki.to_le_bytes())
    .expect("coudn't write to output file");

    // dump current buffer content
    owriter.write_all(data.get_ref()).unwrap();



}

// A number of lines in this function is borrowed from bam2rad funciton at
// https://github.com/COMBINE-lab/alevin-fry/blob/master/libradicl/src/convert.rs
// pub fn read_bamfile(records: &Vec<RadRecord>) {
//     let mut data = Cursor::new(vec![]);

//     for record in records.iter() {
//         tid_list = record.tid_list;
//         if !tid_list.is_empty() {
//             assert!(!tid_list.is_empty(), "Trying to write empty tid_list");
//             let na = tid_list.len();
//             data.write_all(&(na as u32).to_le_bytes()).unwrap();
//             // write bc
//             // data.write_all(&(bc as u32).to_le_bytes()).unwrap();
//             // write umi
//             // data.write_all(&(umi as u32).to_le_bytes()).unwrap();
//             // write tid list
//             for t in tid_list.iter() {
//                 data.write_all(&t.to_le_bytes()).unwrap();
//             }
//         }

//     }
// }