extern crate clap;
extern crate num_cpus;

use clap::{crate_version, App, AppSettings, Arg, ArgGroup};
use std::collections::HashMap;

mod annotation;
mod bam;
mod convert;
mod query_bam_records;
mod rad;

use env_logger::{self, Env};
use log;

fn main() {
    env_logger::Builder::from_env(Env::default().default_filter_or("info")).init();
    let version = crate_version!();
    // let default_num_threads: String = (num_cpus::get() as u32).to_string();
    let default_num_threads = String::from("1");
    let default_max_softlen = String::from("50");
    // let default_supplementary = String::from("keep");
    let app_index = App::new("index")
        .version(version)
        .about("Parse the GTF and build an index to make later runs faster.")
        .arg(Arg::from_usage("-g, --gtf=<FILE> 'Input GTF/GFF file'").required_unless("index").display_order(1))
        .arg(Arg::from_usage("-d, --dir-index=<DIR> 'Output index directory name'").display_order(2))
        .display_order(1);
    let app_bulk = App::new("bulk")
        .version(version)
        .about("Convert alignment of bulk RNA-Seq reads against genome to alignment against transcriptome.")
        .arg(Arg::from_usage("-a, --alignment=<FILE> 'Input SAM/BAM file'").display_order(1))
        .arg(Arg::from_usage("-g, --gtf=<FILE> 'Input GTF/GFF file'").required_unless("index").display_order(1))
        .arg(Arg::from_usage("-i, --index=<DIR> 'Index directory containing parsed GTF files'").required_unless("gtf").display_order(1))
        .arg(Arg::from_usage("-o, --out=<FILE> 'Output file name'").display_order(1))
        .arg(Arg::from_usage("-r, --rad 'Output in RAD format instead of BAM'").display_order(100))
        .arg(Arg::from_usage("-t, --threads=<INT> 'Number of threads for processing bam files'").default_value(&default_num_threads).display_order(2))
        .arg(Arg::from_usage("-s, --max-softclip=<INT> 'Max allowed softclip length'").default_value(&default_max_softlen).display_order(2))
        .group(ArgGroup::with_name("gtf_index_group").args(&["gtf", "index"]).multiple(false).required(true))
        .display_order(2);
        // .arg(Arg::from_usage("--supplementary 'instruction for handling supplementary alignments; one of {keep, keepPrimary, drop}'").default_value(&default_supplementary))
    let app_sc = App::new("sc")
        .version(version)
        .about("Convert alignment of single-cell RNA-Seq reads against genome to alignment against transcriptome.")
        .arg(Arg::from_usage("-a, --alignment=<FILE> 'Input SAM/BAM file'").display_order(1))
        .arg(Arg::from_usage("-g, --gtf=<FILE> 'Input GTF/GFF file'").required_unless("index").display_order(1))
        .arg(Arg::from_usage("-i, --index=<DIR> 'Index directory containing parsed GTF files'").required_unless("gtf").display_order(1))
        .arg(Arg::from_usage("-o, --out=<FILE/DIR> 'Output file name (or directory name when --rad is passed)'").display_order(1))
        .arg(Arg::from_usage("-r, --rad 'Output in RAD format instead of BAM'").display_order(100))
        .arg(Arg::from_usage("-c, --corrected-tags 'Output error-corrected cell barcode and UMI'").display_order(101))
        .arg(Arg::from_usage("-t, --threads=<INT> 'Number of threads for processing bam files'").default_value(&default_num_threads).display_order(2))
        .arg(Arg::from_usage("-s, --max-softclip=<INT> 'Max allowed softclip length'").default_value(&default_max_softlen).display_order(2))
        .arg(Arg::from_usage("-m, --rad-mapped=<FILE> 'Name of output rad file; Only used with --rad'").default_value("map.rad").display_order(3))
        .arg(Arg::from_usage("-u, --rad-unmapped=<FILE> 'Name of file containing the number of unmapped reads for each barcode; Only used with --rad'").default_value("unmapped_bc_count.bin").display_order(3))
        .group(ArgGroup::with_name("gtf_index_group").args(&["gtf", "index"]).multiple(false).required(true))
        .display_order(3);
        // .arg(Arg::from_usage("--supplementary 'instruction for handling supplementary alignments; one of {keep, keepPrimary, drop}'").default_value(&default_supplementary))

    let opts = App::new("mudskipper")
        .setting(AppSettings::SubcommandRequiredElseHelp)
        .setting(AppSettings::DisableHelpSubcommand)
        .version(version)
        .about("Converting RNA-Seq alignments from genome cooridinates to transcriptome coordinates.")
        .subcommand(app_index)
        .subcommand(app_bulk)
        .subcommand(app_sc)
        .get_matches();

    log::info!("Mudskipper started...");
    // convert a SAM/BAM file, in *genome coordinates*,
    // into a BAM file in *transcriptome coordinates*
    if let Some(ref t) = opts.subcommand_matches("index") {
        let ann_file_adr: String = t.value_of("gtf").unwrap().to_string();
        let index_dir: String = t.value_of("dir-index").unwrap().to_string();
        // 
        let mut transcripts_map: HashMap<String, i32> = HashMap::new();
        let mut transcripts: Vec<String> = Vec::new();
        let mut txp_lengths: Vec<i32> = Vec::new();
        annotation::build_tree(&ann_file_adr, &mut transcripts_map, &mut transcripts, &mut txp_lengths, Some(index_dir))
            .expect("cannot build the tree!");
    } else if let Some(ref t) = opts.subcommand_matches("bulk") {
        let bam_file_in: String = t.value_of("alignment").unwrap().to_string();
        let out_file: String = t.value_of("out").unwrap().to_string();
        let threads_count: usize = t.value_of("threads").unwrap().parse::<usize>().unwrap();
        let max_softlen: usize = t.value_of("max-softclip").unwrap().parse::<usize>().unwrap();
        //
        let mut transcripts_map: HashMap<String, i32> = HashMap::new();
        let mut transcripts: Vec<String> = Vec::new();
        let mut txp_lengths: Vec<i32> = Vec::new();
        let trees = if t.is_present("index"){
            let index_dir_adr: String = t.value_of("index").unwrap().to_string();
            annotation::load_tree(&index_dir_adr, &mut transcripts_map, &mut transcripts, &mut txp_lengths)
                .expect("cannot load the tree!")
        } else {
            let ann_file_adr: String = t.value_of("gtf").unwrap().to_string();
            annotation::build_tree(&ann_file_adr, &mut transcripts_map, &mut transcripts, &mut txp_lengths, None)
                .expect("cannot build the tree!")
        };
        if t.is_present("rad") {
            rad::bam2rad_bulk(&bam_file_in, &out_file, &transcripts, &txp_lengths, &trees, &threads_count, &max_softlen);
        } else {
            let required_tags: Vec<&str> = Vec::new();
            bam::bam2bam(
                &bam_file_in,
                &out_file,
                &transcripts,
                &txp_lengths,
                &trees,
                &threads_count,
                &max_softlen,
                &required_tags,
            );
        }
    } else if let Some(ref t) = opts.subcommand_matches("sc") {
        let bam_file_in: String = t.value_of("alignment").unwrap().to_string();
        let out_file: String = t.value_of("out").unwrap().to_string();
        let threads_count: usize = t.value_of("threads").unwrap().parse::<usize>().unwrap();
        let max_softlen: usize = t.value_of("max-softclip").unwrap().parse::<usize>().unwrap();
        let rad_mapped: String = t.value_of("rad-mapped").unwrap().to_string();
        let rad_unmapped: String = t.value_of("rad-unmapped").unwrap().to_string();
        //
        let mut transcripts_map: HashMap<String, i32> = HashMap::new();
        let mut transcripts: Vec<String> = Vec::new();
        let mut txp_lengths: Vec<i32> = Vec::new();
        let trees = if t.is_present("index") {
            let index_dir_adr: String = t.value_of("index").unwrap().to_string();
            annotation::load_tree(&index_dir_adr, &mut transcripts_map, &mut transcripts, &mut txp_lengths)
                .expect("cannot load the tree!")
        } else {
            let ann_file_adr: String = t.value_of("gtf").unwrap().to_string();
            annotation::build_tree(&ann_file_adr, &mut transcripts_map, &mut transcripts, &mut txp_lengths, None)
                .expect("cannot build the tree!")
        };

        let required_tags: Vec<&str>;
        if t.is_present("corrected-tags") {
            required_tags = vec!["CB", "UB"];
        } else {
            required_tags = vec!["CR", "UR"];
        }
        if t.is_present("rad") {
            rad::bam2rad_singlecell(
                &bam_file_in,
                &out_file,
                &rad_mapped,
                &rad_unmapped,
                &transcripts,
                &txp_lengths,
                &trees,
                &threads_count,
                &max_softlen,
                t.is_present("corrected-tags"),
            );
        } else {
            bam::bam2bam(
                &bam_file_in,
                &out_file,
                &transcripts,
                &txp_lengths,
                &trees,
                &threads_count,
                &max_softlen,
                &required_tags,
            );
        }
    }
    log::info!("Mudskipper finished.");
}
