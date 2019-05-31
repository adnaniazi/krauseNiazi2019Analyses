#!/bin/sh
# adjust paths yourself according to where your data is
ALIGNMENT_BAM=/path/to/sorted.aln.bam
GENOME_FA=/path/to/gfp.fa
READS_FA=/path/to/reads.fa
OUTPUT_POLYA_LENGTHS=./path/to/nanopolish_estimates.csv
FAST5_DIR=/path/to/workspace/dir/
SEQUENCING_SUMMARY=/path/to/sequencing_summary.txt
PASS_FASTQ=/path/to/pass_merged.fq

# index the reads
nanopolish index -d $FAST5_DIR -s $SEQUENCING_SUMMARY $PASS_FASTQ

# Convert to fasta
nanopolish extract -r $FAST5_DIR -o $READS_FA

# find poly(A) tails
nanopolish polya --reads $PASS_FASTQ --bam $ALIGNMENT_BAM --genome $GENOME_FA > $OUTPUT_POLYA_LENGTHS
