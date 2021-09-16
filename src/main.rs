extern crate clap;
extern crate num_cpus;

use clap::{crate_version, App, AppSettings, Arg};
// use fnv::FnvHashMap;
// use std::{collections::HashMap, io::BufRead};
use std::collections::HashMap;
use std::env;
// use coitrees::{COITree, IntervalNode};
// use bio_types::strand::Strand;

mod annotation;
mod convert;
mod bam;
mod rad;

// use annotations::ExonNode;

use env_logger;
use log::{info};

fn main() {
    env_logger::init();
    info!("Mudskipper starts...");
    let version = crate_version!();
    let default_num_threads: String = (num_cpus::get() as u32).to_string();
    // let default_num_threads: String = "1".to_string();
    let default_max_softlen: String = "200".to_string();
    let app_bulk = App::new("bulk")
        .version(version)
        .about("Convert alignment of bulk RNA-Seq reads against genome to alignment against transcriptome.")
        .arg(Arg::from("-b, --bam=<bam-file> 'input SAM/BAM file'"))
        .arg(Arg::from("-g, --gtf=<gtf-file> 'input gtf/gff file'"))
        .arg(Arg::from("-o, --out=<output-file> 'output file name'"))
        .arg(Arg::from("-r, --rad 'output in RAD format instead of BAM'"))
        .arg(Arg::from("-t, --threads 'Number of threads for the processing bam files.'").default_value(&default_num_threads))
        .arg(Arg::from("-s, --max-softlen 'Max allowed sofclip length allowed.'").default_value(&default_max_softlen));
    let app_sc = App::new("sc")
        .version(version)
        .about("Convert alignment of single-cell RNA-Seq reads against genome to alignment against transcriptome.")
        .arg(Arg::from("-b, --bam=<bam-file> 'input SAM/BAM file'"))
        .arg(Arg::from("-g, --gtf=<gtf-file> 'input gtf/gff file'"))
        .arg(Arg::from("-o, --out=<output-file> 'output BAM file'"))
        .arg(Arg::from("-t, --threads 'Number of threads for the processing bam files.'").default_value(&default_num_threads))
        .arg(Arg::from("-s, --max-softlen 'Max allowed sofclip length allowed.'").default_value(&default_max_softlen));

    let opts = App::new("mudskipper")
        .setting(AppSettings::SubcommandRequiredElseHelp)
        .setting(AppSettings::DisableHelpSubcommand)
        .version(version)
        .about("Converting RNA-Seq alignments from genome cooridinates to transcriptome coordinates.")
        .subcommand(app_bulk)
        .subcommand(app_sc)
        .get_matches();

    // convert a SAM/BAM file, in *genome coordinates*,
    // into a BAM file in *transcriptome coordinates*
    if let Some(ref t) = opts.subcommand_matches("bulk") {
        let bam_file_in: String = t.value_of_t("bam").unwrap();
        let ann_file_adr: String = t.value_of_t("gtf").unwrap();
        let out_file: String = t.value_of_t("out").unwrap();
        let threads_count: usize = t.value_of_t("threads").unwrap();
        let max_softlen: usize = t.value_of_t("max-softlen").unwrap();
        // 
        let mut transcripts_map: HashMap<String, i32> = HashMap::new();
        let mut transcripts: Vec<String> = Vec::new();
        let mut txp_lengths: Vec<i32> = Vec::new();
        let trees = if std::fs::metadata("parsed_gtf.exon").is_ok() {
            annotation::load_tree(&mut transcripts_map,
                                   &mut transcripts,
                                   &mut &mut txp_lengths).expect("cannot load the tree!")
        } else {
            annotation::build_tree(&ann_file_adr, 
                &mut transcripts_map,
                &mut transcripts,
                &mut txp_lengths).expect("cannot build the tree!")
        };
        if t.is_present("rad") {
            rad::bam2rad(&bam_file_in, &out_file, &transcripts, &txp_lengths, &trees, &threads_count, &max_softlen);
        } else {
            bam::bam2bam(&bam_file_in, &out_file, &transcripts, &txp_lengths, &trees, &threads_count, &max_softlen);
        }
    } else if let Some(ref _t) = opts.subcommand_matches("sc") {
        println!("Support for single-cell bam files not added yet!")
    }
}