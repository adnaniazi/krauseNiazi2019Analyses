library(hdf5r)
library(Biostrings)
library(dplyr)
library(purrr)

demultiplex_reads <- function(read_path=NA,
                              read_id_fast5_file=NA,
                              submat,
                              type,
                              gapOpening,
                              gapExtension,
                              gfp,
                              data='rna',
                              multifast5=F) {
    if (!multifast5){
        f5_obj <- hdf5r::H5File$new(read_path, mode='r')
        f5_tree <- f5_obj$ls(recursive=TRUE)
        f5_tree <- f5_tree$name
        
        raw_read_path <- f5_tree[which(f5_tree == 'Raw/Reads') + 1]
        ri <- f5_obj[[raw_read_path]]$attr_open('read_id')$read()
        bct <- 'Analyses/Basecall_1D_000/BaseCalled_template'
        
        fastq <- f5_obj[[bct]]$open('Fastq')$read()
        #get only the fastq sequence, discarding the quality information and read_id etc
        fastq <-strsplit(fastq, split = "\n")
        fastq <- fastq[[1]][2]
        fastq <- gsub("U", "T", fastq)
        f5_obj$close_all()
        
    } else if (multifast5) {
        filepath <- read_path <- read_id_fast5_file$fast5_file
        f5_obj <- hdf5r::H5File$new(filepath, mode='r')
        read_id <- read_id_fast5_file$read_id
        ri <- read_id
        
        # get eventdata
        event_data_fastq_path <- paste('/', read_id, '/Analyses/Basecall_1D_000/BaseCalled_template/', sep='')
        
        # get the Fastq sequence
        fastq <- f5_obj[[event_data_fastq_path]]$open('Fastq')$read()
        fastq <-strsplit(fastq, split = "\n")
        fastq <- fastq[[1]][2]
        fastq <- gsub("U", "T", fastq)
        f5_obj$close_all()
    }

    ################
    #### RECIPE ####
    ################
    # 1. Check for GFP is found
    # 2. Search for BC 
    # 3. Apply a length filter between 750 and 110
    ################
    
    # Step 1. Check for GFP is found
    num_bases_to_check_gfp_in <- 250
    fastq_start <- substr(fastq, 
                          start=1, 
                          stop=num_bases_to_check_gfp_in)
    fastq_end <- substr(fastq, 
                        start=(nchar(fastq)-num_bases_to_check_gfp_in), 
                        stop=nchar(fastq))
    read_length <- nchar(fastq)
    
    
    as_gfp <- Biostrings::pairwiseAlignment(pattern=gfp,
                                            subject=fastq_start,
                                            substitutionMatrix=submat,
                                            type=type,
                                            scoreOnly=FALSE,
                                            gapOpening=gapOpening,
                                            gapExtension=gapExtension)

    nas_gfp <- as_gfp@score/gfp@length
    
    if (nas_gfp > 0.5 ) {
        read_type <- 'polyA'
    } else {
        read_type <- 'invalid'
        return(list(read_id=ri,
                    filepath=read_path,
                    read_type=read_type,
                    read_length=read_length,
                    read_too_long=NA,
                    read_too_short=NA,
                    nas_gfp=nas_gfp,
                    nas_10bp=NA,
                    nas_30bp=NA,
                    nas_40bp=NA,
                    nas_60bp=NA,
                    nas_100bp=NA,
                    nas_150bp=NA,
                    nas_slip=NA,
                    barcode=NA,
                    barcode_tie=NA,
                    barcode_2=NA,
                    barcode_passed_threshold=NA))
    }
    
    # Step 2: Now search for the barcode
    if (read_type == 'polyA') {
        # search barcode adjacent to the barcode start
        barcode_start <- 1
        barcode_end <- as_gfp@subject@range@start
        
        fastq_to_search_barcode_in <- substr(fastq, start=barcode_start, stop=barcode_end)

        as_all <- Biostrings::pairwiseAlignment(pattern=c('TGCAGATCTCTTGCC', 
                                                          'CTCGAAGCATTGTAA', 
                                                          'AACGGTAGCCACCAA', 
                                                          'TGCACGAGATTGATG', 
                                                          'GACACATAGTCATGG', 
                                                          'CATGAGTGCTGAGCT', 
                                                          'AGCGTGAGACGGCGT'),
                                                subject=fastq_to_search_barcode_in,
                                                substitutionMatrix=submat,
                                                type=type,
                                                scoreOnly=FALSE,
                                                gapOpening=gapOpening,
                                                gapExtension=gapExtension)   
    } 
    
    nas_10bp <- as_all[1]@score/19
    nas_30bp <- as_all[2]@score/19
    nas_40bp <- as_all[3]@score/19
    nas_60bp <- as_all[4]@score/19
    nas_100bp <- as_all[5]@score/19
    nas_150bp <- as_all[6]@score/19
    nas_slip <- as_all[7]@score/19
    
    nas_all <- c(nas_10bp, nas_30bp, nas_40bp, nas_60bp, nas_100bp, nas_150bp, nas_slip)
    bc_threshold <- 0.6
    
    m <- which(nas_all==max(nas_all))
    barcode_tie <- F
    barcode_2 <- NA
    barcode_passed_threshold <- F
    count <- 1
    for (i in m) {
        if (nas_all[i] != 0) {
            if (count == 1) {
                count = count + 1
                if (nas_all[i] > bc_threshold) barcode_passed_threshold = T
                if (i == 1) barcode <- 10
                else if (i == 2) barcode <- 30
                else if (i == 3) barcode <- 40
                else if (i == 4) barcode <- 60
                else if (i == 5) barcode <- 100
                else if (i == 6) barcode <- 150
                else if (i == 7) barcode <- 999
            } else if (count == 2) {
                count <- count + 1
                barcode_tie <- T
                if (i == 1) barcode_2 <- 10
                else if (i == 2) barcode_2 <- 30
                else if (i == 3) barcode_2 <- 40
                else if (i == 4) barcode_2 <- 60
                else if (i == 5) barcode_2 <- 100
                else if (i == 6) barcode_2 <- 150
                else if (i == 7) barcode_2 <- 999
            } else break  
        } else {
            barcode <- 0
        }
        
    }
    
    # 3. Apply a length filter between 750 and 1100
    rtl <- FALSE
    rts <- FALSE
    if (nchar(fastq) > 1100 ) {
        rtl <- TRUE
    } else if (nchar(fastq) < 700 ) {
        rts <- TRUE
    }
    
    return(list(read_id=ri,
                filepath=read_path,
                read_type=read_type,
                read_length=read_length,
                read_too_long=rtl,
                read_too_short=rts,
                nas_gfp=nas_gfp,
                nas_10bp=nas_10bp,
                nas_30bp=nas_30bp,
                nas_40bp=nas_40bp,
                nas_60bp=nas_60bp,
                nas_100bp=nas_100bp,
                nas_150bp=nas_150bp,
                nas_slip=nas_slip,
                barcode=barcode,
                barcode_tie=barcode_tie,
                barcode_2=barcode_2,
                barcode_passed_threshold=barcode_passed_threshold))
}


demultiplex_reads_parallel <- function(fast5_dir, num_cores, save_path, multifast5=F, data='rna'){
    # Define the data
    # make gfp sequence
    gfp <- Biostrings::DNAString('CCACCATGGTGAGCAAGGGCGAGGAGCTG')
    rc_gfp <- Biostrings::reverseComplement(gfp)
    # Define a scoring matrix
    match <- 1
    mismatch <- -1
    type <- 'local'
    gapOpening <- 0
    gapExtension <- 1
    submat <- Biostrings::nucleotideSubstitutionMatrix(match=match,
                                                       mismatch=mismatch,
                                                       baseOnly=TRUE)

    # search for all the fast5 files in the user-specified directory
    message('\t- Searching for all Fast5 files...\r')
    fast5_files_list <- list.files(path=fast5_dir,
                                   pattern="\\.fast5$",
                                   recursive=TRUE,
                                   full.names=TRUE)
    num_files <- length(fast5_files_list)
    message('\t  Done! Found ', num_files, ' fast5 files\r')
    
    message('\t- Starting a parallel cluster...\r')
    # Initiate cluster
    cl <- parallel::makeCluster(num_cores) #, outfile='')
    doSNOW::registerDoSNOW(cl)
    `%dopar%` <- foreach::`%dopar%`
    `%do%` <- foreach::`%do%`
    
    on.exit(parallel::stopCluster(cl))
    
    message('\t  Done!')
    
    if (!multifast5) {    
        # Split the data into chunks
        files_per_chunk <- 4000
        total_files <- length(fast5_files_list)
        total_chunks <- ceiling(total_files/files_per_chunk)
        
        counter <- 0
        result <- list()
        #loop
        message('\t- Searching for barcodes...\r')
        for(chunk in c(1:total_chunks)){
            if(chunk == total_chunks)
                fast5_files_subset <- fast5_files_list[((counter*files_per_chunk)+1):total_files]
            else
                fast5_files_subset <- fast5_files_list[((counter*files_per_chunk)+1):((counter+1)*files_per_chunk)]
            counter <- counter + 1
            
            # progress bar
            message('\r') 
            
            message(paste('\t  Processing chunk ', chunk, ' of ', total_chunks, '\r', sep='')) 
            pb <- txtProgressBar(min=1, max=length(fast5_files_subset), style=3)
            progress <- function(n) setTxtProgressBar(pb, n)
            opts <- list(progress=progress)
            
            mcoptions <- list(preschedule=TRUE, set.seed=FALSE, cleanup=TRUE)
            data_list <-foreach::foreach(file_path = fast5_files_subset,
                                         .combine = 'rbind',
                                         .options.snow = opts,
                                         .inorder=FALSE,
                                         .export='demultiplex_reads',
                                         .packages=c('Biostrings', 'hdf5r'),
                                         .options.multicore = mcoptions) %dopar% {
                                             tryCatch({
                                                 demultiplex_reads(read_path=file_path,
                                                                   read_id_fast5_file=NA,
                                                                   submat=submat,
                                                                   type=type,
                                                                   gapOpening=gapOpening,
                                                                   gapExtension=gapExtension,
                                                                   gfp=gfp,
                                                                   data=data,
                                                                   multifast5=multifast5)
                                             },
                                             error=function(e){
                                                 print(e)
                                                 ls <- list(read_id=NA,
                                                            filepath=file_path,
                                                            read_type=NA,
                                                            read_length=NA,
                                                            read_too_long=NA,
                                                            read_too_short=NA,
                                                            nas_gfp=NA,
                                                            nas_10bp=NA,
                                                            nas_30bp=NA,
                                                            nas_40bp=NA,
                                                            nas_60bp=NA,
                                                            nas_100bp=NA,
                                                            nas_150bp=NA,
                                                            nas_slip=NA,
                                                            barcode=NA,
                                                            barcode_tie=NA,
                                                            barcode_2=NA,
                                                            barcode_passed_threshold=NA)
                                             })
                                         }
            result[[chunk]] <- data_list
        }
    } else if (multifast5) {
        message('\t- Discovering reads in the ',  num_files, ' multifast5 files found...\r')
        read_id_fast5_file <- dplyr::tibble(read_id=character(), fast5_file=character())
        for (fast5_file in fast5_files_list) {
            f5_obj <- hdf5r::H5File$new(fast5_file, mode='r')
            f5_tree <- f5_obj$ls(recursive=F)
            f5_tree <- f5_tree$name
            f5_tree <- dplyr::mutate(dplyr::tbl_df(f5_tree), fast5_file=fast5_file)
            f5_tree <- dplyr::rename(f5_tree, read_id=value)
            read_id_fast5_file <- rbind(read_id_fast5_file, f5_tree)
            f5_obj$close_all()
        }
        message('\t  Done! Found ', nrow(read_id_fast5_file), ' reads\r')
        # convert the data frame to list with rows as elements of the list
        read_id_fast5_file <- split(read_id_fast5_file, seq(nrow(read_id_fast5_file)))
        
        # Split the data into chunks
        files_per_chunk <- 4000
        total_files <- length(read_id_fast5_file)
        total_chunks <- ceiling(total_files/files_per_chunk)

        #loop
        message('\t- Searching for barcodes ...\r')
        counter <- 0
        result <- list()
        for(chunk in c(1:total_chunks)){
            # divide data in chunks
            if(chunk == total_chunks)
                read_id_fast5_file_subset <- read_id_fast5_file[((counter*files_per_chunk)+1):total_files]
            else
                read_id_fast5_file_subset <- read_id_fast5_file[((counter*files_per_chunk)+1):((counter+1)*files_per_chunk)]
            counter <- counter + 1
            message(paste('\t  Processing chunk ', chunk, ' of ', total_chunks, '\r', sep=''))
            
            # progress bar
            message('\r')
            pb <- txtProgressBar(min=1, max=length(read_id_fast5_file_subset), style=3)
            progress <- function(n) setTxtProgressBar(pb, n)
            opts <- list(progress=progress)
            # foreach loop
            mcoptions <- list(preschedule=TRUE, set.seed=FALSE, cleanup=TRUE)
            data_list <- foreach::foreach(riff=read_id_fast5_file_subset,
                                          .combine='rbind',
                                          .inorder=FALSE,
                                          .errorhandling='pass',
                                          .options.snow=opts,
                                          .export='demultiplex_reads',
                                          .packages=c('Biostrings', 'hdf5r'),
                                          .options.multicore=mcoptions) %dopar% {
                                              tryCatch({
                                                  demultiplex_reads(read_path=file_path,
                                                                    read_id_fast5_file=riff,
                                                                    submat=submat,
                                                                    type=type,
                                                                    gapOpening=gapOpening,
                                                                    gapExtension=gapExtension,
                                                                    gfp=gfp,
                                                                    data=data,
                                                                    multifast5=multifast5)
                                              },
                                              error=function(e){
                                                  print(e)
                                                  ls <- list(read_id=NA,
                                                             filepath=file_path,
                                                             read_type=NA,
                                                             read_length=NA,
                                                             read_too_long=NA,
                                                             read_too_short=NA,
                                                             nas_gfp=NA,
                                                             nas_10bp=NA,
                                                             nas_30bp=NA,
                                                             nas_40bp=NA,
                                                             nas_60bp=NA,
                                                             nas_100bp=NA,
                                                             nas_150bp=NA,
                                                             nas_slip=NA,
                                                             barcode=NA,
                                                             barcode_tie=NA,
                                                             barcode_2=NA,
                                                             barcode_passed_threshold=NA)
                                              })
                                          }
            result[[chunk]] <- data_list
        }
    }
    
    message('\r')
    # Close the cluster
    message('\t- Shutting down the cluster...\r')
    close(pb)

    # Unlist and make a dataframe
    message('\t- Making a dataframe.\r')
    
    result <- purrr::map(result, function(.x) tibble::as_tibble(.x))
    result <- dplyr::bind_rows(result, .id = "chunk")
    result <- dplyr::select(result, -chunk)
    message('\t- Saving the table in CSV file.\r')
    data.table::fwrite(result, save_path)
    message('\t  All done! Now quitting.\r')
    
    return(result)
}

# please adjust the paths/parameters as per your own setup
num_cores <- 120
data <- 'rna'
multifast5 <- F

fast5_dir <- '/export/valenfs/data/processed_data/MinION/20181217_max_rna_polyashield_spike/basecalled_data/workspace/pass'
save_path <- '/export/valenfs/data/processed_data/MinION/20181217_max_rna_polyashield_spike/polya_estimation/decoded_barcode.csv'
df <- demultiplex_reads_parallel(fast5_dir, num_cores, save_path, multifast5=multifast5, data=data)




