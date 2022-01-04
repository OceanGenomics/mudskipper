# mudskipper

## Quick start
```bash
# convert a bulk RNA-Seq genomic BAM to a transcriptomic BAM for quantification with Salmon
mudkipper bulk --gtf annotation.gtf --alignment genomic.bam --out transcriptomic.bam
# convert a single-cell RNA-Seq genomic SAM to a transcriptomic RAD for quantification with alevin-fry
mudkipper sc --gtf annotation.gtf --alignment genomic.sam --out transcriptomic_dir

# build and store a GTF index; useful if you want to convert multiple BAM/SAM files
mudskipper index --gtf annotation.gtf --dir-index gtf_index
# run mudskipper using the pre-built GTF index
mudkipper bulk --index gtf_index --alignment genomic.bam --out transcriptomic.bam
mudkipper sc --index gtf_index --alignment genomic.sam --out transcriptomic_dir
```

## Introduction

`mudskipper` is a tool for converting genomic BAM/SAM files to transcriptomic BAM/RAD files. More specifically, it projects the genomic coordinates of each alignment entry to transcriptomic coordinates based on a given transcript annotation in GTF format.

The current focus in `mudskipper` is to enable transcript quantification from genomic alignments without re-mapping short reads onto the transcriptome (e.g. using `salmon` for bulk RNA-Seq samples and `alevin-fry` for single-cell RNA-Seq samples). The motivation for developing `mudskipper` was that there are many tools that require alignment of short RNA-Seq reads against reference genome. Transcript quantification tools on the other hand, usually expect alignment of short RNA-Seq reads against the transcriptome. `mudskipper` enables users to perform transcript quantification using genomic alignments instead of starting the whole process from scratch (by avoiding mapping short RNA-Seq reads directly against the transcriptome).

**Algorithmic details:** `mudskipper` parses the input GTF file and builds an interval tree that stores the coordinates of exons for *all* the transcripts present in the GTF file. It then processes each genomic alignment of the given input BAM/SAM file. Note that genomic alignments of RNA-Seq reads are spliced (i.e. they might have `N` in the CIGAR string). Therefore, `mudskipper` uses the generated interval tree to project each genomic alignment to all transcripts that satisfy the following conditions:
1. Transcripts whose start and end coordinates on the genome mark a region that fully contains the alignment,
2. Transcripts whose exon coordinates match the splicing of the genomic alignment

For each such transcript, a BAM record is created that stores the proper alignment to that transcript. **Note** that the first condition above means that projection to intergenic or intronic regions are not reported. It also means that `mudskipper` currently does not report support reads that align transcripts in an overhanging fashion. This requirement can be relaxed in the future to support overhanging alignments.

## Building from source

`mudskipper` is written in [Rust](https://www.rust-lang.org/) and can be built as follows:
[**Requires:** [Rust toolchain](https://www.rust-lang.org/tools/install)]
```bash
cargo build --release
```

This will create a binary executable at `target/release/mudskipper`. You can add this file to your environment PATH variable for convenience. This can be done temporarily using:
```bash
export PATH=`pwd`/target/release/:$PATH
```

## Interface

### Projection of bulk RNA-Seq read alignments
```
mudskipper bulk [OPTIONS] --alignment <FILE> --out <FILE> (--gtf <FILE>|--index <DIR>)

OPTIONS:
    -a, --alignment <FILE>     Input SAM/BAM file
    -g, --gtf <FILE>           Input GTF/GFF file
    -i, --index <DIR>          Index directory containing parsed GTF files
    -o, --out <FILE>           Output file name
    -s, --max-softlen <INT>    Max allowed softclip length [default: 200]
    -t, --threads <INT>        Number of threads for processing bam files [default: 1]
    -r, --rad                  Output in RAD format instead of BAM
    -h, --help                 Prints help information
    -V, --version              Prints version information
```

### Projection of single-cell RNA-Seq read alignments
```
mudskipper sc [OPTIONS] --alignment <FILE> --out <FILE/DIR> (--gtf <FILE>|--index <DIR>)

OPTIONS:
    -a, --alignment <FILE>       Input SAM/BAM file
    -g, --gtf <FILE>             Input GTF/GFF file
    -i, --index <DIR>            Index directory containing parsed GTF files
    -o, --out <FILE/DIR>         Output file name (or directory name when --rad is passed)
    -s, --max-softlen <INT>      Max allowed softclip length [default: 200]
    -t, --threads <INT>          Number of threads for processing bam files [default: 1]
    -m, --rad-mapped <FILE>      The name of output rad file; Only used with --rad [default: map.rad]
    -u, --rad-unmapped <FILE>    The name of file containing the number of unmapped reads for each barcode; Only used with --rad [default: unmapped_bc_count.bin]
    -r, --rad                    Output in RAD format instead of BAM
    -c, --corrected-tags         Output error-corrected cell barcode and UMI
    -h, --help                   Prints help information
    -V, --version                Prints version information
```

### Build and store the GTF interval tree
```
mudskipper index --dir-index <DIR> --gtf <FILE>

OPTIONS:
    -g, --gtf <FILE>         Input GTF/GFF file
    -d, --dir-index <DIR>    Output index directory name
    -h, --help               Prints help information
    -V, --version            Prints version information
```
