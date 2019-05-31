library(tailfindr)
fast5_dir <- '/specify/path/'
save_dir <- '/specify/path/'
num_cores <- 120
csv_filename <- 'rna_tails.csv'
df <- find_tails(fast5_dir=fast5_dir,
                 save_dir=save_dir,
                 csv_filename=csv_filename,
                 num_cores=num_cores)



