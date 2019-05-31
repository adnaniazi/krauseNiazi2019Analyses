library(here)
# source functions
source(here('extract-read-data.R'))
source(here('get-moves-per-read.R'))
source(here('get-total-moves-between-borders.R'))

data_path <- here('data') # adjust accordinly to where the data is

### Step 1: get read paths for ff and ab, and tail start and ends for both, and read_id
dna_ff_path <- file.path(data_path, 'dna-krause-lsk108_109-flipflop_basecalling-tailfindr_estimates-two_replicates-with_filepaths.csv')
dna_st_path <- file.path(data_path, 'dna-krause-lsk108_109-standard_basecalling-tailfindr_estimates-two_replicates-with_filepaths.csv')

dna_ff_df <- read.csv(dna_ff_path, header = TRUE, stringsAsFactors = FALSE)
dna_st_df <- read.csv(dna_st_path, header = TRUE, stringsAsFactors = FALSE)

dna_ff_paths <- dna_ff_df$file_path
dna_st_paths <- dna_st_df$file_path

tail_start_ff <- dna_ff_df$tail_start
tail_end_ff <- dna_ff_df$tail_end
read_id_ff <- dna_ff_df$read_id

tail_start_st <- dna_st_df$tail_start
tail_end_st <- dna_st_df$tail_end
read_id_st <- dna_st_df$read_id


# process flipflop data
# Ebarassingly slow code; no time to parallalize it
data_list_ff <- list()
d = 0
for (i in seq(length(dna_ff_paths))) {
    cat('Flipflop read #: ', i, '\n')
    if (is.na(tail_start_ff[i]) | is.na(tail_end_ff[i])) {
        d <- list(read_id = read_id_ff[i],
                  moves_in_tail = NA)
        data_list_ff[[i]] <- d
        next()
    }
    
    tryCatch({
        d <- get_moves_per_read(file_path = dna_ff_paths[i],
                                basecalled_with = 'guppy',
                                multifast5 = FALSE,
                                model = 'flipflop',
                                tail_start = tail_start_ff[i],
                                tail_end = tail_end_ff[i])
        },
        error=function(e){
            d <- list(read_id = read_id_ff[i],
                      moves_in_tail = NA)
        })
    data_list_ff[[i]] <- d
}

# process albacore data
# Ebarassingly slow code; no time to parallalize it

data_list_st <- list()
for (i in seq(length(dna_st_paths))) {
    cat('Albacore read #: ', i, '\n')
    if (is.na(tail_start_st[i]) | is.na(tail_end_st[i])) {
        d <- list(read_id = read_id_st[i],
                  moves_in_tail = NA)
        data_list_st[[i]] <- d
        next()
    }
    
    tryCatch({
        d <- get_moves_per_read(file_path = dna_st_paths[i],
                                basecalled_with = 'albacore',
                                multifast5 = FALSE,
                                model = 'standard',
                                tail_start = tail_start_st[i],
                                tail_end = tail_end_st[i])
    },
    error=function(e){
        d <- list(read_id = read_id_st[i],
                  moves_in_tail = NA)
    })
    data_list_st[[i]] <- d
}

### Step 3: combine all informtion into a data frame and save it
st_moves_df <- purrr::map(data_list_st, function(.x) tibble::as_tibble(.x))
st_moves_df <- dplyr::bind_rows(st_moves_df, .id = "chunk")
st_moves_df <- dplyr::select(st_moves_df, -chunk)

ff_moves_df <- purrr::map(data_list_ff, function(.x) tibble::as_tibble(.x))
ff_moves_df <- dplyr::bind_rows(ff_moves_df, .id = "chunk")
ff_moves_df <- dplyr::select(ff_moves_df, -chunk)

# renames ab ff columns
st_moves_df <- dplyr::rename(st_moves_df, moves_in_tail_st = moves_in_tail)
ff_moves_df <- dplyr::rename(ff_moves_df, moves_in_tail_ff = moves_in_tail)

# save the dataframe
data.table::fwrite(st_moves_df, file.path(data_path, 'dna-krause-lsk108_109-standard_basecalling-moves_in_tail-two_replicates-with_filepaths.csv'))
data.table::fwrite(ff_moves_df, file.path(data_path, 'dna-krause-lsk108_109-flipflop_basecalling-moves_in_tail-two_replicates-with_filepaths.csv'))
