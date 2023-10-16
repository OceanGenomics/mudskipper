#!/bin/bash
set -ex

cut -f 1,3 Aligned.toTranscriptome.out.sam | grep -v '^@' | sort | uniq > read.transcript.STAR.txt

cut -f 1,3 Aligned.mudskipperout.sam | grep -v '^@' | sort | uniq > read.transcript.mudskipper.txt  

join -v 1 read.transcript.STAR.txt read.transcript.mudskipper.txt > only_star.txt  