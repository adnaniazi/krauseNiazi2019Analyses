library(tailfindr)
df_dna <- find_tails(fast5_dir = /specify/fast5_dir,
                     save_dir = /specify/fast5_dir,
                     basecall_group = 'Basecall_1D_000',
                     csv_filename = 'tails.csv',
                     num_cores = 120,
                     save_plots = FALSE,
                     dna_datatype = 'pcr-dna')


