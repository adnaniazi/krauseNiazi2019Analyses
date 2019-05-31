#!/bin/sh
# please adjust the path according to your own setup
source /Home/ii/adnann/nanopore/teloprime/env36-kjmepefuru/bin/activate

INPUT=/export/valenfs/data/raw_data/minion/20190515_DNA_Max_polyA-DNA-spikes2/single_fast5/
OUTPUT=/export/valenfs/data/processed_data/MinION/20190515_DNA_Max_polyA-DNA-spikes2/basecalled_data_albacore

nice -n 0 read_fast5_basecaller.py \
--flowcell FLO-MIN106 \
--kit SQK-LSK108 \
--recursive \
--output_format fast5,fastq \
--input $INPUT \
--save_path $OUTPUT \
--worker_threads 12 \
2> logfile.log
