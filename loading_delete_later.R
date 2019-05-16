drake::loadd(rna001_kr_st_bc)
drake::loadd(rna001_kr_st_tf)
drake::loadd(rna001_kr_st_np)
drake::loadd(rna001_kr_st_tr)
decoded_barcodes = rna001_kr_st_bc
tailfindr_estimates = rna001_kr_st_tf
nanopolish_estimates = rna001_kr_st_np
transcript_start_info = rna001_kr_st_tr


drake::loadd(rna002_kr_st_bc)
drake::loadd(rna002_kr_st_tf)
drake::loadd(rna002_kr_st_np)
drake::loadd(rna002_kr_st_tr)
decoded_barcodes = rna002_kr_st_bc
tailfindr_estimates = rna002_kr_st_tf
nanopolish_estimates = rna002_kr_st_np
transcript_start_info = rna002_kr_st_tr


rm(list = ls())
drake::loadd(dna_kr_ff_bc)
drake::loadd(dna_kr_st_tf)
drake::loadd(dna_kr_ff_tf)
drake::loadd(dna_kr_st_tr)
decoded_barcodes_ff = dna_kr_ff_bc
tailfindr_estimates_st = dna_kr_st_tf
tailfindr_estimates_ff = dna_kr_ff_tf
transcript_start_info_st = dna_kr_st_tr
