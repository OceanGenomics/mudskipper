extern crate clap;

use clap::{crate_authors, crate_version, App, AppSettings, Arg, ArgSettings};
use std::collections::HashMap;
use std::env;

mod annotations;
mod intersection;
mod bam;

fn main() {
    let crate_authors = crate_authors!("\n");
    let version = crate_version!();

    let bam2bam_app = App::new("bam2bam")
        .version(version)
        .author(crate_authors)
        .about("Convert genome aligned bam/sam file to transcriptome aligned bam file.")
        .arg(Arg::from("-b, --bam=<bam-file> 'input SAM/BAM file'"))
        .arg(Arg::from("-g, --gtf=<gtf-file> 'input gtf/gff file'"))
        .arg(Arg::from("-o, --out=<output-file> 'output BAM file'"));

    let opts = App::new("mudskipper")
        .setting(AppSettings::SubcommandRequiredElseHelp)
        .version(version)
        .author(crate_authors)
        .about("Converting the genome aligned cooridinates to transcriptome coordinates.")
        .subcommand(bam2bam_app)
        .get_matches();

    // convert a SAM/BAM file, in *genome coordinates*,
    // into a BAM file in *transcriptome coordinates*
    if let Some(ref t) = opts.subcommand_matches("bam2bam") {
        let bam_file_in: String = t.value_of_t("bam").unwrap();
        let ann_file_adr: String = t.value_of_t("gtf").unwrap();
        let bam_file_out: String = t.value_of_t("out").unwrap();
        //libradicl::convert::bam2bam(input_file, rad_file, num_threads, &log)

        //let ann_file_adr: String = env::args().nth(1).expect("missing src");
        let mut transcripts_map: HashMap<String, i32> = HashMap::new();
        let mut transcripts: Vec<String> = Vec::new();
        let mut txp_lengths: Vec<i32> = Vec::new();
        let trees = annotations::build_tree(&ann_file_adr, 
                                            &mut transcripts_map,
                                            &mut transcripts,
                                            &mut txp_lengths).expect("cannot build the tree!");
        //let bam_file_in: String = env::args().nth(2).expect("missing src");
        //let bam_file_out: String = env::args().nth(3).expect("missing src");
        bam::read_bamfile(&bam_file_in, &bam_file_out, &transcripts, &txp_lengths, &trees);
    }
}