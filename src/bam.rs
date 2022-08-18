use coitrees::COITree;

use rust_htslib::bam::{header, Format, Header, Writer};

use crate::annotation;
use annotation::ExonNode;

extern crate bio_types;

use crate::convert;
use crate::query_bam_records::BAMQueryRecordReader;

extern crate fnv;
use fnv::FnvHashMap;

#[allow(clippy::too_many_arguments)]
pub fn bam2bam(
    input_bam_filename: &String,
    output_bam_filename: &String,
    transcripts: &[String],
    txp_lengths: &[i32],
    trees: &FnvHashMap<String, COITree<ExonNode, u32>>,
    threads_count: &usize,
    max_softlen: &usize,
    max_overhang: &usize,
    required_tags: &[&str],
) -> i32 {
    // setup the output BAM Writer
    let mut output_header = Header::new();
    // output_header.push_record(header::HeaderRecord::new(b"HD").push_tag(b"VN", &"1.4"));
    log::info!("Number of reference sequences: {}", transcripts.len());
    for (counter, tname) in transcripts.iter().enumerate() {
        let tlen = txp_lengths[counter];
        output_header.push_record(header::HeaderRecord::new(b"SQ").push_tag(b"SN", &tname).push_tag(b"LN", &tlen));
    }
    let mut output_writer = Writer::from_path(output_bam_filename, &output_header, Format::Bam).unwrap();

    let reader_threads: Option<usize>;
    if *threads_count >= 2 {
        let threads_count_half = threads_count / 2;
        log::info!("thread count: {} in total", threads_count_half * 2);
        log::info!("thread count: {} for reading", threads_count_half);
        log::info!("thread count: {} for writing", threads_count_half);
        reader_threads = Some(threads_count_half);
        output_writer
            .set_threads(threads_count_half)
            .expect("Failed to set number of BAM writing threads.");
    } else {
        reader_threads = None;
        log::info!("thread count: {}", threads_count);
    }

    // setup the input BAM Reader
    let mut bqr = BAMQueryRecordReader::new(input_bam_filename, reader_threads);
    let input_header = bqr.get_header().to_owned();

    let mut missed = 0i32;

    loop {
        match bqr.get_next_query_records() {
            Ok(ret_vec) => match ret_vec {
                Some(ret_vec) => {
                    if ret_vec.is_empty() {
                        missed += 1;
                    }
                    for r in ret_vec.iter() {
                        let txp_records = convert::convert_query_bam_records(
                            r,
                            &input_header,
                            transcripts,
                            txp_lengths,
                            trees,
                            max_softlen,
                            max_overhang,
                            required_tags,
                        );
                        if txp_records.is_empty() {
                            missed += 1;
                        }
                        for txp_rec in txp_records.iter() {
                            output_writer.write(txp_rec).unwrap();
                        }
                    }
                }
                None => break,
            },
            Err(e) => {
                log::error!(
                    "Cannot find the mate for query {}. The mate might be missing or not in the right order! \
                If the sam/bam file is not name sorted, please consider using \"mudskipper shuffle\" or adding the flag \"-l\". \
                If you wish to skip the unpaired query, please consider adding the flag \"-p\".",
                    e
                );
                panic!("Invalid input sam/bam file!")
            }
        }
    }
    missed
}

// skip the unpaired query for paired-end reads
#[allow(clippy::too_many_arguments)]
pub fn bam2bam_skip(
    input_bam_filename: &String,
    output_bam_filename: &String,
    transcripts: &[String],
    txp_lengths: &[i32],
    trees: &FnvHashMap<String, COITree<ExonNode, u32>>,
    threads_count: &usize,
    max_softlen: &usize,
    max_overhang: &usize,
    required_tags: &[&str],
) -> i32 {
    // setup the output BAM Writer
    let mut output_header = Header::new();
    // output_header.push_record(header::HeaderRecord::new(b"HD").push_tag(b"VN", &"1.4"));
    log::info!("Number of reference sequences: {}", transcripts.len());
    for (counter, tname) in transcripts.iter().enumerate() {
        let tlen = txp_lengths[counter];
        output_header.push_record(header::HeaderRecord::new(b"SQ").push_tag(b"SN", &tname).push_tag(b"LN", &tlen));
    }
    let mut output_writer = Writer::from_path(output_bam_filename, &output_header, Format::Bam).unwrap();

    let reader_threads: Option<usize>;
    if *threads_count >= 2 {
        let threads_count_half = threads_count / 2;
        log::info!("thread count: {} in total", threads_count_half * 2);
        log::info!("thread count: {} for reading", threads_count_half);
        log::info!("thread count: {} for writing", threads_count_half);
        reader_threads = Some(threads_count_half);
        output_writer
            .set_threads(threads_count_half)
            .expect("Failed to set number of BAM writing threads.");
    } else {
        reader_threads = None;
        log::info!("thread count: {}", threads_count);
    }

    // setup the input BAM Reader
    let mut bqr = BAMQueryRecordReader::new(input_bam_filename, reader_threads);
    let input_header = bqr.get_header().to_owned();

    let mut missed = 0i32;

    while let Some(ret_vec) = bqr.get_next_query_records_skip() {
        if ret_vec.is_empty() {
            missed += 1;
        }
        for r in ret_vec.iter() {
            let txp_records = convert::convert_query_bam_records(
                r,
                &input_header,
                transcripts,
                txp_lengths,
                trees,
                max_softlen,
                max_overhang,
                required_tags,
            );
            if txp_records.is_empty() {
                missed += 1;
            }

            for txp_rec in txp_records.iter() {
                output_writer.write(txp_rec).unwrap();
            }
        }
    }
    missed
}

// NOTE: this was an attempt at re-producing this bug: https://github.com/OceanGenomics/mudskipper/issues/6
//       But failed.
// pub fn test_read_bam() -> () {
//     let mut input_bam = Reader::from_path(String::from("test.bam")).unwrap();
//     input_bam.set_threads(2).expect("Failed to set number of reading threads.");
//     // for r in input_bam.records() {
//     //     let rec = r.expect("Failed parsing the BAM/SAM file");
//     //     println!("qname:{} tid:{} tpos:{}", String::from_utf8(rec.qname().to_vec()).unwrap(), rec.tid(), rec.pos());
//     // }
//     let mut record = Record::new();
//     while let Some(result) = input_bam.read(&mut record) {
//         result.expect("Failed to parse BAM record");
//         println!("qname:{} tid:{} tpos:{}", String::from_utf8(record.qname().to_vec()).unwrap(), record.tid(), record.pos());
//     }
// }
