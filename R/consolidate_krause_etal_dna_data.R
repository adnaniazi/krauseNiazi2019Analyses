#' Title
#'
#' @param decoded_barcodes_ff 
#' @param tailfindr_estimates_ff 
#' @param transcript_start_info_ff 
#' @param tailfindr_estimates_st 
#' @param transcript_start_info_st 
#' @importFrom magrittr %>% 
#' @importFrom magrittr %<>% 
#' @return
#' @export
#'
#' @examples
consolidate_krause_etal_dna_data <- function(decoded_barcodes_ff,
                                             tailfindr_estimates_ff,
                                             transcript_start_info_ff,
                                             tailfindr_estimates_st,
                                             transcript_start_info_st) {
  
  decoded_barcodes_ff %<>% 
    # keep only valid reads in which GFP sequence is there
    dplyr::filter(read_type == 'polyA' | read_type == 'polyT') %>% 
    # remove reads that are shorter than 700bp or greater than 900bp
    # such reads have a high chance that the barcode is either not 
    # present or that the read is fused with any other read -- making 
    # the barcode -- if detecded -- unreliable
    dplyr::filter(read_length >= 700 & read_length <= 900) %>% 
    # remove reads in which two or more barcodes have an equally
    # high alignment score. We don't need such ambiguous reads in our
    # analysis
    dplyr::filter(!barcode_tie) %>% 
    # reomve reads in which the barcode didn't pass
    dplyr::filter(barcode_passed_threshold) %>%  
    # remove slip oligo (barcode == 999)
    dplyr::filter(barcode != 999) %>% 
    # select only read_id and barcode column as we will need only
    # these two columns subsequently
    dplyr::select(read_id, read_type, barcode, replicate, file_path)
    
  tailfindr_estimates_ff %<>% 
    # filter out reads for which tailfindr's tail length
    # estimates are not available
    dplyr::filter(!is.na(tail_length)) %>% 
    # take only valid reads
    dplyr::filter(tail_is_valid) %>% 
    # remove irrelevant columns
    dplyr::select(-read_type, -tail_is_valid, -file_path, -replicate) %>% 
    # rename variables:
    dplyr::rename(tail_start_ff = tail_start,
                  tail_end_ff = tail_end,
                  samples_per_nt_ff = samples_per_nt,
                  tail_length_ff = tail_length)
  
  tailfindr_estimates_st %<>% 
    # filter out reads for which tailfindr's tail length
    # estimates are not available
    dplyr::filter(!is.na(tail_length)) %>% 
    # take only valid reads
    dplyr::filter(tail_is_valid == TRUE) %>% 
    # remove irrelevant columns
    dplyr::select(-read_type, -tail_is_valid, -file_path, -replicate) %>% 
    # rename variables:
    dplyr::rename(tail_start_st = tail_start,
                  tail_end_st = tail_end,
                  samples_per_nt_st = samples_per_nt,
                  tail_length_st = tail_length)
    
  transcript_start_info_st %<>% 
    dplyr::select(-file_path) %>% 
    dplyr::rename(transcript_alignment_start_st = transcript_alignment_start)
  
  transcript_start_info_ff %<>% 
    dplyr::select(-file_path) %>% 
    dplyr::rename(transcript_alignment_start_ff = transcript_alignment_start)
  
  # join all these four datasets
  df <- dplyr::left_join(tailfindr_estimates_ff, tailfindr_estimates_st, by = "read_id") %>% 
    dplyr::left_join(decoded_barcodes_ff, by = "read_id") %>% 
    dplyr::left_join(transcript_start_info_st, by = "read_id") %>% 
    dplyr::left_join(transcript_start_info_ff, by = "read_id") 
  
  # remove rows with missing tail length predictions
  df <- df[!is.na(df$tail_length_ff), ]
  df <- df[!is.na(df$tail_length_st), ]
  # remove rows with missing barcode predictions
  df <- df[!is.na(df$barcode), ]
  df <- dplyr::distinct(df)
  
  # make factor columns for categoical variables
  df %<>% dplyr::mutate(read_type = as.factor(read_type)) %>% 
    dplyr::mutate(barcode = as.factor(barcode)) %>% 
    dplyr::mutate(replicate = as.factor(replicate)) 
  
  return(df)
}
