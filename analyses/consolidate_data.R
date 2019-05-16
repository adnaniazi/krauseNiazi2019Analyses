rna001_kr_data = consolidate_krause_etal_rna001_data(
  decoded_barcodes = rna001_kr_st_bc,
  tailfindr_estimates = rna001_kr_st_tf,
  nanopolish_estimates = rna001_kr_st_np,
  transcript_start_info = rna001_kr_st_tr
)

rna002_kr_data = consolidate_krause_etal_rna002_data(
  decoded_barcodes = rna002_kr_st_bc,
  tailfindr_estimates = rna002_kr_st_tf,
  nanopolish_estimates = rna002_kr_st_np,
  transcript_start_info = rna002_kr_st_tr
)

rna_kr_data = combine_rna001_rna002_datasets(
  rna001_dataset = rna001_kr_data,
  rna002_dataset = rna002_kr_data
)

rna_wo_data = consolidate_workman_etal_data(
  tailfindr_estimates = rna_wo_st_tf,
  nanopolish_estimates = rna_wo_st_np,
  transcript_start_info = rna_wo_st_tr
  ) 
  
dna_kr_data = consolidate_krause_etal_dna_data(
  decoded_barcodes_ff = dna_kr_ff_bc,
  tailfindr_estimates_ff = dna_kr_ff_tf,
  transcript_start_info_ff = dna_kr_ff_tr,
  tailfindr_estimates_st = dna_kr_st_tf,
  transcript_start_info_st = dna_kr_st_tr
)