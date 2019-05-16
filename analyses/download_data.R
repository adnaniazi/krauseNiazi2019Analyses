# This script downloads all the data into extdata directory

# ------------- Description of abrreviations used -------------------
# bc = barcode, kr = krause et al data,  ff = flipflop, 
# bc = decoded barcodes, tf = tailfindr predictions,
# np = nanopolish predictions, wo = workman et al data
# dna = DNA data, rna001 = RNA data obtained using RNA001 kit
# rna002 = RNA data obtained using RNA002 kit
# -------------------------------------------------------------------

#####################################################################
# ------------ Download Krause et al DNA data ----------------------#
# Data contains two replicates: one on LSK108 kit and the other on
# LSL109 kit. Apart from containing spikeins of known lengths, the 
# data also contains reads from zebrafish transcriptome collected as
# part of the experiment. These reads will be filtered out in
# downstream processing steps
#####################################################################

# get barcode information based on data basecalled using flipflop model
dna_kr_ff_bc =
  get_data(release_version = 'v1.0.2',
           file_name = 'dna-krause-lsk108_109-flipflop_basecalling-decoded_barcodes-two_replicates-with_filepaths.csv')

# get tailfindr predictions on data basecalled using flipflop model
dna_kr_ff_tf =
  get_data(release_version = 'v1.0.2',
           file_name = 'dna-krause-lsk108_109-flipflop_basecalling-tailfindr_estimates-two_replicates-with_filepaths.csv')

# get tailfindr predictions on data basecalled using standrard model
dna_kr_st_tf =
  get_data(release_version = 'v1.0.2',
           file_name = 'dna-krause-lsk108_109-standard_basecalling-tailfindr_estimates-two_replicates-with_filepaths.csv')

# get transcript_start using alignment with eGFP on data basecalled using standrard model
dna_kr_st_tr =
  get_data(release_version = 'v1.0.2',
           file_name = 'dna-krause-lsk108_109-standard_basecalling-transcript_alignment_start-two_replicates-with_filepaths.csv')

# get transcript_start using alignment with eGFP on data basecalled using flipflop model
dna_kr_ff_tr =
  get_data(release_version = 'v1.0.2',
           file_name = 'dna-krause-lsk108_109-flipflop_basecalling-transcript_alignment_start-two_replicates-with_filepaths.csv')

#####################################################################
# ------------ Download Krause et al RNA001 data -------------------#
# Data contains two replicates: both obtained using RNA001 kit. 
# Apart from containing spikeins of known lengths, the 
# data also contains reads from zebrafish transcriptome collected as
# part of the experiment. These reads will be filtered out in
# downstream processing steps
#####################################################################

# get barcode information based on data basecalled using standard model
rna001_kr_st_bc =
  get_data(release_version = 'v1.0.2',
           file_name = 'rna-krause-rna001-standard_basecalling-decoded_barcodes-two_replicates-with_filepaths.csv')

# get nanopolish predictions on data basecalled using standrard model
rna001_kr_st_np =
  get_data(release_version = 'v1.0.2',
           file_name = 'rna-krause-rna001-standard_basecalling-nanopolish_estimates-two_replicates-with_filepaths.csv')

# get tailfindr predictions on data basecalled using standrard model
rna001_kr_st_tf =
  get_data(release_version = 'v1.0.2',
           file_name = 'rna-krause-rna001-standard_basecalling-tailfindr_estimates-two_replicates-with_filepaths.csv')

# get transcript_start using alignment with eGFP
rna001_kr_st_tr =
  get_data(release_version = 'v1.0.2',
           file_name = 'rna-krause-rna001-standard_basecalling-transcript_alignment_start-two_replicates-with_filepaths.csv')

#####################################################################
# ------------------ Download RNA002 data --------------------------#
# Contain data obtained using RNA002 kit. 
#####################################################################

# get barcode information based on data basecalled using standard model
rna002_kr_st_bc =
  get_data(release_version = 'v1.0.2',
           file_name = 'rna-krause-rna002-standard_basecalling-decoded_barcodes-with_filepaths.csv')

# get nanopolish predictions on data basecalled using standrard model
rna002_kr_st_np =
  get_data(release_version = 'v1.0.2',
           file_name = 'rna-krause-rna002-standard_basecalling-nanopolish_estimates-with_filepaths.csv')

# get tailfindr predictions on data basecalled using standrard model
rna002_kr_st_tf =
  get_data(release_version = 'v1.0.2',
           file_name = 'rna-krause-rna002-standard_basecalling-tailfindr_estimates-with_filepaths.csv')

# get transcript_start using alignment with eGFP
rna002_kr_st_tr =
  get_data(release_version = 'v1.0.2',
           file_name = 'rna-krause-rna002-standard_basecalling-transcript_alignment_start-with_filepaths.csv')

# get barcode information based on data basecalled using fliflop model
rna002_kr_ff_bc =
  get_data(release_version = 'v1.0.2',
           file_name = 'rna-krause-rna002-flipflop_basecalling-decoded_barcodes-with_filepaths.csv')

# get nanopolish predictions on data basecalled using fliflop model
rna002_kr_ff_np =
  get_data(release_version = 'v1.0.2',
           file_name = 'rna-krause-rna002-flipflop_basecalling-nanopolish_estimates-with_filepaths.csv')

# get tailfindr predictions on data basecalled using fliflop model
rna002_kr_ff_tf =
  get_data(release_version = 'v1.0.2',
           file_name = 'rna-krause-rna002-flipflop_basecalling-tailfindr_estimates-with_filepaths.csv')

#####################################################################
# ------------ Download Workman et al RNA001 data ------------------#
# The data contains polyA tails of known lengths. The data was 
# generated using RNA001 kit.
#####################################################################

# get nanopolish predictions on data basecalled using standrard model
rna_wo_st_np =
  get_data(release_version = 'v1.0.2',
           file_name = 'rna-workman-rna001-standard_basecalling-nanopolish_estimates-all_datasets-with_filepaths.csv')

# get tailfindr predictions on data basecalled using standrard model
rna_wo_st_tf =
  get_data(release_version = 'v1.0.2',
           file_name = 'rna-workman-rna001-standard_basecalling-tailfindr_estimates-all_datasets-with_filepaths.csv')

# get transcript_start using alignment with enolase
rna_wo_st_tr =
  get_data(release_version = 'v1.0.2',
           file_name = 'rna-workman-rna001-standard_basecalling-transcript_alignment_start-all_datasets-with_filepaths.csv')
