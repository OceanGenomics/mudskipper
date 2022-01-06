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

## Table of Contents
- [Quick start](#quick-start)
- [Introduction](#introduction)
- [Building from source](#building-from-source)
- [Interface](#interface)
  - [Projection of bulk RNA-Seq read alignments](#projection-of-bulk-rna-seq-read-alignments)
  - [Projection of single-cell RNA-Seq read alignments](#projection-of-single-cell-rna-seq-read-alignments)
  - [Building and storing the GTF interval tree](#building-and-storing-the-gtf-interval-tree)
- [Limitations](#limitations)

## Introduction

`mudskipper` is a tool for converting genomic BAM/SAM files to transcriptomic BAM/RAD files. More specifically, it projects the genomic coordinates of each alignment entry to transcriptomic coordinates based on a given transcript annotation in GTF format.

The current focus in `mudskipper` is to enable transcript quantification from genomic alignments without re-mapping short reads onto the transcriptome (e.g. using [`salmon`](https://github.com/COMBINE-lab/salmon) for bulk RNA-Seq samples and [`alevin-fry`](https://github.com/COMBINE-lab/alevin-fry/) for single-cell RNA-Seq samples). The motivation for developing `mudskipper` was that there are many tools that require alignment of short RNA-Seq reads against reference genome. Transcript quantification tools on the other hand, usually expect alignment of short RNA-Seq reads against the transcriptome. `mudskipper` enables users to perform transcript quantification using genomic alignments instead of starting the whole process from scratch (by avoiding mapping short RNA-Seq reads directly against the transcriptome).

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
mudskipper bulk [OPTIONS] --alignment <FILE> (--gtf <FILE>|--index <DIR>) --out <FILE>

OPTIONS:
    -a, --alignment <FILE>      Input SAM/BAM file
    -g, --gtf <FILE>            Input GTF/GFF file
    -i, --index <DIR>           Index directory containing parsed GTF files
    -o, --out <FILE>            Output file name
    -s, --max-softclip <INT>    Max allowed softclip length [default: 50]
    -t, --threads <INT>         Number of threads for processing bam files [default: 1]
    -r, --rad                   Output in RAD format instead of BAM
    -h, --help                  Prints help information
    -V, --version               Prints version information
```

#### Required arguments
##### `-a, --alignment <FILE>`
Using this argument you can specify the input alignment file. Currenlty BAM/SAM formats are supported. This BAM/SAM file should contain alignment of short RNA-Seq reads against reference genome.
> ✏️ The alignments stored in this file are potentially spliced. 

##### `-g, --gtf <FILE>`
This argument can be used to pass the gene annotation in GTF format. In this case, the interval tree is build from the GTF file on the fly. Alternatively, [`--index`](#-i---index-dir) argument can be used. That means `--gtf` and `--index` are mutually exclusive.
> ✏️ Make sure that the GTF file corresponds to the same version of the reference genome to which short reads have been aligned. If some target sequences are missing from the GTF file, alignments to those target sequences will be dropped automatically.

##### `-i, --index <DIR>`
This argument specifies the path to the pre-built interval tree, previously created from a GTF file (see [Building and storing the GTF interval tree](#building-and-storing-the-gtf-interval-tree) for more details about how to create such an index). This is helpful if you wish to run `mudskipper` on many BAM/SAM files.  the  Note that `--gtf` and `--index` are mutually exclusive.

##### `-o, --out <FILE>`
The path to the output alignment file. By default, the output alignment file is in BAM format. If [`--rad`](#-r---rad) is passed, then the output alignment file will be in RAD format.

#### Optional arguments
##### `-r, --rad`
Pass this argument to output alignments in RAD format instead of BAM. This argument is not set by default.

##### `-t, --threads <INT>`
The number of threads to use for reading and writing BAM files. By default, this is set to 1.

##### `-s, --max-softclip <INT>`
Drop alignments that have more than INT soft-clipped bases. By default, this is set to 50.

### Projection of single-cell RNA-Seq read alignments
```
mudskipper sc [OPTIONS] --alignment <FILE> (--gtf <FILE>|--index <DIR>) --out <FILE/DIR>

OPTIONS:
    -a, --alignment <FILE>       Input SAM/BAM file
    -g, --gtf <FILE>             Input GTF/GFF file
    -i, --index <DIR>            Index directory containing parsed GTF files
    -o, --out <FILE/DIR>         Output file name (or directory name when --rad is passed)
    -s, --max-softclip <INT>     Max allowed softclip length [default: 50]
    -t, --threads <INT>          Number of threads for processing bam files [default: 1]
    -r, --rad                    Output in RAD format instead of BAM
    -m, --rad-mapped <FILE>      The name of output rad file; Only used with --rad [default: map.rad]
    -u, --rad-unmapped <FILE>    The name of file containing the number of unmapped reads for each barcode; Only used with --rad [default: unmapped_bc_count.bin]
    -c, --corrected-tags         Output error-corrected cell barcode and UMI
    -h, --help                   Prints help information
    -V, --version                Prints version information
```

#### Required arguments
##### `-a, --alignment <FILE>`
Using this argument you can specify the input alignment file. Currenlty BAM/SAM formats are supported. This BAM/SAM file should contain alignment of short RNA-Seq reads against reference genome.
> ✏️ The alignments stored in this file are potentially spliced. 

##### `-g, --gtf <FILE>`
This argument can be used to pass the gene annotation in GTF format. In this case, the interval tree is build from the GTF file on the fly. Alternatively, [`--index`](#-i---index-dir-1) argument can be used. That means `--gtf` and `--index` are mutually exclusive.
> ✏️ Make sure that the GTF file corresponds to the same version of the reference genome to which short reads have been aligned. If some target sequences are missing from the GTF file, alignments to those target sequences will be dropped automatically.

##### `-i, --index <DIR>`
This argument specifies the path to the pre-built interval tree, previously created from a GTF file (see [Building and storing the GTF interval tree](#building-and-storing-the-gtf-interval-tree) for more details about how to create such an index). This is helpful if you wish to run `mudskipper` on many BAM/SAM files.  the  Note that `--gtf` and `--index` are mutually exclusive.

##### `-o, --out <FILE/DIR>`
The path to the output alignment file in BAM format. If [`--rad`](#-r---rad-1) is passed, this argument specifies the output directory that contains the RAD format as well as other files required by `alevin-fry` for performing transcript quantification.

#### Optional arguments
##### `-r, --rad`
Pass this argument to output alignments in RAD format instead of BAM. This argument is not set by default.

##### `-m, --rad-mapped <FILE>`
Specifies the name of output rad file. This option is used only when [`--rad`](#-r---rad-1) is passed. In that case, this file will be stored in the directory passed by [`--out`](#-o---out-filedir). By default, this is set to `map.rad`.

##### `-u, --rad-unmapped <FILE>`
Specifies the name of the file containing the number of unmapped reads for each barcode. This option is used only when [`--rad`](#-r---rad-1) is passed. In that case, this file will be stored in the directory passed by [`--out`](#-o---out-filedir). By default, this is set to `unmapped_bc_count.bin`.

##### `-t, --threads <INT>`
The number of threads to use for reading and writing BAM files. By default, this is set to 1.

##### `-s, --max-softclip <INT>`
Drop alignments that have more than INT soft-clipped bases. By default, this is set to 50.

##### `-c, --corrected-tags`
By default, `mudskipper sc` expects to see `CR` and `UR` tags in the input BAM/SAM file, which store cell barcode and UMI respectively. If this argument is passed, then `CB` and `UB` tags must be present instead, which store corrected cell barcode and corrected UMI respectively.

### Building and storing the GTF interval tree
```
mudskipper index --gtf <FILE> --dir-index <DIR>

OPTIONS:
    -g, --gtf <FILE>         Input GTF/GFF file
    -d, --dir-index <DIR>    Output index directory name
    -h, --help               Prints help information
    -V, --version            Prints version information
```

When using `mudskipper bulk` or `mudskipper sc`, if you pass the gene annotation file, `mudskipper` builds an interval tree on the fly. This interval tree is used to query transcriptomic coordinates that overlap a given genomic region. Instead of building the interval tree on the fly, you can build the interval tree and store it for later use. This is helpful if you wish to run `mudskipper` on many BAM/SAM files.

#### Required arguments
##### `-g, --gtf <FILE>`
Specifies the gene/transcript annotation file in GTF format.

##### `-d, --dir-index <DIR>`
The path of the directory where the interval tree files will be stored.
> ✏️ This directory will be created if does not exist.

## Limitations
`mudskipper` is still in early stages of development with lots of room for improvements. So far, `mudskipper` has been tested only for the purpose of transcript quantification. Currently, it has the following known limitations:
- Chimeric alignments are not reported. That means for now a projected alignment is reported only if all its segments fall on the same target.
- It only reports projected alignments that are fully contained in a transcript. In other words, it currenlty does not report any overhanging alignments. [[#10](/../../issues/10)]
- Some fields and optional tags of the output BAM might not be properly updated. [[#13](/../../issues/13)]
- For single-cell samples, it drops alignment of reads that have `N` in their barcode. [[#15](/../../issues/15)]

Bug reports and suggestions are warmly welcomed. 
