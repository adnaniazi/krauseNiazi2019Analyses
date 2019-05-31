#!/bin/sh

# please adjust the path as per your own setup
source /Home/ii/adnann/nanopore/teloprime/env36-kjmepefuru/bin/activate

INPUT=/export/valenfs/data/raw_data/minion/20181217_max_RNA_polyAshield-spike/fast5
OUTPUT=/export/valenfs/data/processed_data/MinION/20181217_max_rna_polyashield_spike/basecalled_data


nice -n 0 read_fast5_basecaller.py \
	--flowcell FLO-MIN106 \
	--kit SQK-RNA001 \
	--recursive \
	--output_format fast5,fastq \
	--input $INPUT \
	--save_path $OUTPUT \
	--worker_threads 12 \
	2> logfile.log
