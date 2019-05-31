#!/bin/sh

# please adjust the path acoording to your own setup
INPUT=/export/valenfs/data/raw_data/minion/20190515_DNA_Max_polyA-DNA-spikes2/single_fast5
OUTPUT=/export/valenfs/data/processed_data/MinION/20190515_DNA_Max_polyA-DNA-spikes2/basecalled_data_flipflop

/Home/ii/adnann/kjempefuru_software/guppy_231/ont-guppy-cpu/bin/./guppy_basecaller \
-c dna_r9.4.1_450bps_flipflop.cfg \
-i $INPUT \
-s $OUTPUT \
-r \
--fast5_out \
--hp_correct 1 \
--disable_pings 1 \
--enable_trimming 0 \
--num_callers 120 \


