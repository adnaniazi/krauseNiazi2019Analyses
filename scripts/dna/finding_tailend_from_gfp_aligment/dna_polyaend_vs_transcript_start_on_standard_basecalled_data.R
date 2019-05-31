library(hdf5r)
library(Biostrings)
library(dplyr)
library(purrr)

find_sample_index_for_fastq_base <-
    function(event_data, fastq_base_number) {
        model_state_length <- nchar(event_data$model_state[1])
        cumm_sum <- cumsum(event_data$move) + model_state_length - 1
        event_data <- cbind(event_data, cumm_sum)
        
        row_index <- min(which(cumm_sum == fastq_base_number))
        if (row_index == -Inf | row_index == Inf) {
            row_index <- min(which(cumm_sum == (fastq_base_number - 1)))
        }
        
        sample_index <- event_data$start[row_index]
        return(sample_index)
    }


find_transcript_start <- function(file_path = NA,
                                  read_id_fast5_file = NA,
                                  read_type = 'invalid',
                                  basecall_group = 'Basecall_1D_000',
                                  submat,
                                  type,
                                  gapOpening,
                                  gapExtension,
                                  gfp_end,
                                  data = 'rna',
                                  multifast5 = FALSE) {
    if (!multifast5) {
        f5_obj <- hdf5r::H5File$new(file_path, mode='r')
        f5_tree <- f5_obj$ls(recursive=TRUE)
        f5_tree <- f5_tree$name
        # define all the paths
        raw_signal_path <- grep('.*Signal$', f5_tree, perl = TRUE, value = TRUE)
        if (sum(grepl('Raw/Reads/Read_[0-9]+$', f5_tree)) == 0) {
            read_id_path <-  grep('.*Raw$', f5_tree, perl = TRUE, value = TRUE)
        } else {
            read_id_path <- f5_tree[grepl('Raw/Reads/Read_[0-9]+$', f5_tree)]
        }
        # make fastq, Event/Move path
        event_data_fastq_path <- grep(paste0('.*', basecall_group, '/BaseCalled_template$'),
                                      f5_tree, perl = TRUE, value = TRUE)
        # make segmentation path based on the basecall group
        sp <- strsplit(basecall_group, split = '_')
        seg_group <- sp[[1]][3]
        segmentation_path <- grep(paste0('.*', seg_group, '/Summary/segmentation$'),
                                  f5_tree, perl = TRUE, value = TRUE)
        # make basecalled_template path
        basecall_1d_template_path <- grep(paste0('.*', basecall_group, '/Summary/basecall_1d_template$'),
                                          f5_tree, perl = TRUE, value = TRUE)
    }
    else {
        f5_obj <- hdf5r::H5File$new(read_id_fast5_file$fast5_file, mode='r')
        full_read_id <- read_id_fast5_file$read_id
        # define all the paths
        raw_signal_path <- paste('/', full_read_id, '/Raw/Signal', sep='')
        event_data_fastq_path <- paste0('/', full_read_id, '/Analyses/', basecall_group, '/BaseCalled_template/')
        read_id_path <- paste('/', full_read_id, '/Raw', sep='')

        # make segmentation path based on the basecall group
        sp <- strsplit(basecall_group, split='_')
        seg_group <- sp[[1]][3]
        segmentation_path <- paste0('/', full_read_id, '/Analyses/Segmentation_', seg_group, '/Summary/segmentation')
        basecall_1d_template_path <- paste0('/', full_read_id, '/Analyses/', basecall_group, '/Summary/basecall_1d_template')
    }

    # get the data
    raw_data <- f5_obj[[raw_signal_path]]$read()
    read_id <- f5_obj[[read_id_path]]$attr_open('read_id')$read()
    fastq <- f5_obj[[event_data_fastq_path]]$open('Fastq')$read()
    start <- f5_obj[[segmentation_path]]$attr_open('first_sample_template')$read()
    called_events <- f5_obj[[basecall_1d_template_path]]$attr_open('called_events')$read()
    # compute called_event if it has been set to W by ONT in the FAST5 file (Wierd! I know)
    compute_called_events <- ifelse(called_events == 'W', TRUE, FALSE)
    # get the event data, if present -- or make it, if not present
    make_event_data = FALSE
    if ('Events' %in% names(f5_obj[[event_data_fastq_path]])){
        event_data <- f5_obj[[event_data_fastq_path]]$open('Events')$read()
        if (!('length' %in% colnames(event_data))) {
            stride <- f5_obj[[basecall_1d_template_path]]$attr_open('block_stride')$read()
        } else {
            stride <- event_data$length[1]
        }
    } else {
        move <- f5_obj[[event_data_fastq_path]]$open('Move')$read()
        stride <- f5_obj[[basecall_1d_template_path]]$attr_open('block_stride')$read()
        make_event_data = TRUE
    }

    # extract just the fastq removing quality scores
    fastq <- strsplit(fastq, split = "\n")
    fastq <- fastq[[1]][2]


    # if event_data wasn't present, make it now
    if (make_event_data) {
        # The line below is there to remove R CMD CHECK
        # "no visible binding for global variable" error
        move_cumsum <- fastq_bases <- NULL

        start_col <- seq(from = start,
                         to = start + stride * (length(move) - 1),
                         by = stride)
        event_data <- data.frame(move = move,
                                 start = start_col,
                                 move_cumsum = cumsum(move),
                                 fastq_bases = fastq,
                                 stringsAsFactors = FALSE)

        # The line below is there to remove R CMD CHECK
        # "no visible binding for global variable" error
        model_state <- NULL
        event_data <- dplyr::mutate(event_data,
                                    model_state = substr(fastq_bases,
                                                         start=move_cumsum,
                                                         stop=move_cumsum))
        event_data <- dplyr::select(event_data, model_state, start, move)
    }
    
    # if event_data wasn't present, make it now
    if (make_event_data) {
        # The line below is there to remove R CMD CHECK
        # "no visible binding for global variable" error
        move_cumsum <- fastq_bases <- NULL
        
        start_col <- seq(
            from = start,
            to = start + stride * (length(move) - 1),
            by = stride
        )
        event_data <- data.frame(
            move = move,
            start = start_col,
            move_cumsum = cumsum(move),
            fastq_bases = fastq,
            stringsAsFactors = FALSE
        )
        
        # The line below is there to remove R CMD CHECK
        # "no visible binding for global variable" error
        model_state <- NULL
        event_data <- dplyr::mutate(
            event_data,
            model_state = substr(fastq_bases,
                                 start = move_cumsum,
                                 stop = move_cumsum)
        )
        event_data <-
            dplyr::select(event_data, model_state, start, move)
    }
    
    transcript_start = NA
    
    if (read_type == 'polyA') {
      #############################
      #### RECIPE POLY-A Reads ####
      #############################
      # 1. Reverse both the fastq and the gfp end to check
      # 2. Align the reversed gfp with the revered fastq
      # 3. If the alignment starts with index 1, then you have found the perfect
      #    candidate read to compare the transcript start site with the poly(A)
      #    end site
      ############################
      
      # 1. Reverse both the fastq and the gfp end to check
      fastq_rev <- Biostrings::reverse(fastq)
      gfp_end_rev <- Biostrings::reverse(gfp_end)
      
      # 2. Align reverse of gfp with reverse of fastq
      num_bases_to_check_gfp_in <- 75
      fastq_segment_to_check_gfp_in <- substr(fastq_rev,
                                              start = 1,
                                              stop = num_bases_to_check_gfp_in)
      
      as_gfp <- Biostrings::pairwiseAlignment(
        pattern = gfp_end_rev,
        subject = fastq_segment_to_check_gfp_in,
        substitutionMatrix = submat,
        type = type,
        scoreOnly = FALSE,
        gapOpening = gapOpening,
        gapExtension = gapExtension
      )
      
      nas_gfp <- as_gfp@score / gfp_end_rev@length
      
      # Step 2: Now search for the barcode
      transcript_start <- NA
      if (nas_gfp > 0.5) {
        alignment_gfp_start_index <- as_gfp@pattern@range@start
        alignment_fastq_start_index <- nchar(fastq) - as_gfp@subject@range@start + 2
        
        if (alignment_gfp_start_index == 1 &
            as.character(substr(as_gfp@pattern, 1, 3)) == as.character(substr(gfp_end_rev, 1 , 3)) & 
            as.character(substr(as_gfp@subject, 1, 3)) == as.character(substr(gfp_end_rev, 1 , 3))) {
          sample_index <-
            find_sample_index_for_fastq_base (
              event_data = event_data,
              fastq_base_number = alignment_fastq_start_index
            )
          transcript_start <- sample_index
        }
      }
    } else if (read_type == 'polyT') {
      #############################
      #### RECIPE POLY-T Reads ####
      #############################
      # 1. Reverse complement GFP end
      # 2. Align the reversed gfp with the revered fastq
      # 3. If the alignment starts with index 1, then you have found the perfect
      #    candidate read to compare the transcript start site with the poly(A)
      #    end site
      ############################
      
      # 1. Reverse complement GFP end
      gfp_end_rev_comp <- Biostrings::reverseComplement(gfp_end)
      
      # 2. Align reverse complement of gfp with the fastq
      num_bases_to_check_gfp_in <- 75
      fastq_segment_to_check_gfp_in <- substr(fastq,
                                              start = 1,
                                              stop = num_bases_to_check_gfp_in)
      
      as_gfp <- Biostrings::pairwiseAlignment(
        pattern = gfp_end_rev_comp,
        subject = fastq_segment_to_check_gfp_in,
        substitutionMatrix = submat,
        type = type,
        scoreOnly = FALSE,
        gapOpening = gapOpening,
        gapExtension = gapExtension
      )
      
      nas_gfp <- as_gfp@score / gfp_end_rev_comp@length
      
      # Step 2: Now search for the barcode
      transcript_start <- NA
      if (nas_gfp > 0.5) {
        alignment_gfp_start_index <- as_gfp@pattern@range@start
        alignment_fastq_start_index <- as_gfp@subject@range@start + 1
        
        if (alignment_gfp_start_index == 1 &
            as.character(substr(as_gfp@pattern, 1, 3)) == as.character(substr(gfp_end_rev_comp, 1 , 3)) &
            as.character(substr(as_gfp@subject, 1, 3)) == as.character(substr(gfp_end_rev_comp, 1 , 3))) {
          sample_index <-
            find_sample_index_for_fastq_base (
              event_data = event_data,
              fastq_base_number = alignment_fastq_start_index
            )
          transcript_start <- sample_index
        }
      }
    }
    return(list(transcript_start = transcript_start,
                read_id = read_id,
                file_path = file_path))
}

# # Define the data
# # make gfp sequence
# gfp_end <-Biostrings::DNAString('ATCACTCTCGGCATGGACGAGCTGTACAAGTAG')
# # Define a scoring matrix
# match <- 1
# mismatch <- -1
# type <- 'local'
# gapOpening <- 0
# gapExtension <- 1
# submat <- Biostrings::nucleotideSubstitutionMatrix(match = match,
#                                                    mismatch = mismatch,
#                                                    baseOnly = TRUE)
# polyt_read_1 <- "/export/valenfs/data/processed_data/MinION/20181219_max_dna_pcrspikes/basecalled_data/workspace/pass/36/sars_HP_Z240_Tower_Workstation_20181219_FAK36778_MN21607_sequencing_run_PCR_48535_read_80896_ch_365_strand.fast5"
# polyt_read_2 <- "/export/valenfs/data/processed_data/MinION/20181219_max_dna_pcrspikes/basecalled_data/workspace/pass/36/sars_HP_Z240_Tower_Workstation_20181219_FAK36778_MN21607_sequencing_run_PCR_48535_read_81530_ch_365_strand.fast5"
# 
# polya_read_1 <- "/export/valenfs/data/processed_data/MinION/20181219_max_dna_pcrspikes/basecalled_data/workspace/pass/35/sars_HP_Z240_Tower_Workstation_20181219_FAK36778_MN21607_sequencing_run_PCR_48535_read_43283_ch_308_strand.fast5"
# polya_read_2 <- "/export/valenfs/data/processed_data/MinION/20181219_max_dna_pcrspikes/basecalled_data/workspace/pass/35/sars_HP_Z240_Tower_Workstation_20181219_FAK36778_MN21607_sequencing_run_PCR_48535_read_59932_ch_69_strand.fast5"
# 
# df <- find_transcript_start(
#   file_path = polyt_read_2,
#   read_id_fast5_file = NA,
#   read_type = 'polyT',
#   basecall_group = 'Basecall_1D_000',
#   submat,
#   type,
#   gapOpening,
#   gapExtension,
#   gfp_end,
#   data = 'dna',
#   multifast5 = FALSE
# )
# print(df)


find_transcript_start_parallel <-
    function(fast5_files_list,
             read_types_list,
             num_cores,
             multifast5 = F,
             data = 'rna') {
        # Define the data
        # make gfp sequence
        gfp_end <-
            Biostrings::DNAString('ATCACTCTCGGCATGGACGAGCTGTACAAGTAG')
        # Define a scoring matrix
        match <- 1
        mismatch <- -1
        type <- 'local'
        gapOpening <- 0
        gapExtension <- 1
        submat <-
            Biostrings::nucleotideSubstitutionMatrix(match = match,
                                                     mismatch = mismatch,
                                                     baseOnly = TRUE)

        # search for all the fast5 files in the user-specified directory
        message('\t- Searching for all Fast5 files...\r')
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
            total_chunks <- ceiling(total_files / files_per_chunk)

            counter <- 0
            result <- list()
            #loop
            message('\t- Searching for barcodes...\r')
            for (chunk in c(1:total_chunks)) {
                if (chunk == total_chunks) {
                    fast5_files_subset <-
                      fast5_files_list[((counter * files_per_chunk) + 1):total_files]
                    read_types_subset <-
                      read_types_list[((counter * files_per_chunk) + 1):total_files]
                } else {
                    fast5_files_subset <-
                      fast5_files_list[((counter * files_per_chunk) + 1):((counter + 1) * files_per_chunk)]
                    read_types_subset <-
                      read_types_list[((counter * files_per_chunk) + 1):((counter + 1) * files_per_chunk)]
                }
                counter <- counter + 1

                # progress bar
                message('\r')

                message(paste('\t  Processing chunk ', chunk, ' of ', total_chunks, '\r', sep = ''))
                
                pb <-
                    txtProgressBar(
                        min = 1,
                        max = length(fast5_files_subset),
                        style = 3
                    )
                progress <- function(n)
                    setTxtProgressBar(pb, n)
                opts <- list(progress = progress)

                mcoptions <-
                    list(
                        preschedule = TRUE,
                        set.seed = FALSE,
                        cleanup = TRUE
                    )
                data_list <-
                    foreach::foreach(
                        file_path = fast5_files_subset,
                        read_type = read_types_subset,
                        .combine = 'rbind',
                        .options.snow = opts,
                        .inorder = FALSE,
                        .export = c(
                            'find_transcript_start',
                            'find_sample_index_for_fastq_base'
                        ),
                        .packages = c('Biostrings', 'hdf5r'),
                        .options.multicore = mcoptions
                    ) %dopar% {
                        tryCatch({
                            find_transcript_start(
                                file_path = file_path,
                                read_id_fast5_file = NA,
                                read_type = read_type,
                                basecall_group = 'Basecall_1D_000',
                                submat = submat,
                                type = type,
                                gapOpening = gapOpening,
                                gapExtension = gapExtension,
                                gfp_end = gfp_end,
                                data = 'rna',
                                multifast5 = FALSE
                            )
                        },
                        error = function(e) {
                            print(e)
                            ls <-
                                list(
                                    transcript_start = NA,
                                    read_id = NA,
                                    file_path = file_path
                                )
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

        result <- purrr::map(result, function(.x)
            tibble::as_tibble(.x))
        result <- dplyr::bind_rows(result, .id = "chunk")
        result <- dplyr::select(result,-chunk)
        message('\t  All done! Now quitting.\r')

        return(result)
    }

## Please adjust the path here accordingly
base_path <- '/valenfs/piggyback_data_for_upload/data_with_filepaths'
dna_data_path <-
    file.path(base_path, 'dna-krause-lsk108_109-standard_basecalling-tailfindr_estimates-two_replicates-with_filepaths.csv')
dna_data <- read.csv(dna_data_path, header = TRUE, stringsAsFactors = FALSE)
num_cores <- 120
data <- 'dna'
multifast5 <- FALSE

fast5_files_list <- dna_data$file_path
read_types_list <- dna_data$read_type

df <- find_transcript_start_parallel(fast5_files_list = fast5_files_list,
                                     read_types_list = read_types_list,
                                     num_cores = num_cores,
                                     multifast5 = multifast5,
                                     data = data)

df <- df %>% dplyr::rename(transcript_alignment_start = transcript_start)
df <- dplyr::bind_rows(df, .id = "chunk")
df <- dplyr::select(df, -chunk)
# cleanup the tibble
df <- tidyr::unnest(df)

save_path <- file.path(base_path, 'dna-krause-lsk108_109-standard_basecalling-transcript_alignment_start-two_replicates-with_filepaths.csv')
message('\t- Saving the table in CSV file.\r')
data.table::fwrite(df, save_path)
