#' Title
#'
#' @param tailfindr_estimates 
#' @param nanopolish_estimates 
#' @param transcript_start_info 
#'
#' @importFrom magrittr %>% 
#'
#' @return
#' @export
#'
#' @examples
consolidate_workman_etal_data <- function(tailfindr_estimates,
                                          nanopolish_estimates,
                                          transcript_start_info) {

  tailfindr_estimates <- tailfindr_estimates %>% 
    # filter out reads for which tailfindr's tail length
    # estimates are not available
    dplyr::filter(!is.na(tail_length)) %>% 
    # remove irrelevant columns
    dplyr::select(-file_path) %>% 
    # rename variables:
    dplyr::rename(tail_start_tf = tail_start,
                  tail_end_tf = tail_end,
                  samples_per_nt_tf = samples_per_nt,
                  tail_length_tf = tail_length)
    
  nanopolish_estimates <- nanopolish_estimates %>% 
    # take only PASS reads
    dplyr::filter(qc_tag == 'PASS') %>% 
    # reomve irrelevant columns
    dplyr::select(-contig, -position, -leader_start, 
                  -adapter_start, -file_path, -dataset) %>% 
    # calculate samples_per_nt rate from nanopolish's read rate
    dplyr::mutate(samples_per_nt_np = 3012/read_rate) %>% 
    # rename variables for consistency sake
    dplyr::rename(read_id = readname,
                  tail_start_np = polya_start,
                  tail_end_np = transcript_start,
                  tail_length_np = polya_length,
                  read_rate_np = read_rate,
                  qc_tag_np = qc_tag)
    
  transcript_start_info <- transcript_start_info  %>% 
    dplyr::select(-file_path)
  
  df <- dplyr::left_join(tailfindr_estimates, nanopolish_estimates, by = "read_id") %>% 
    dplyr::left_join(transcript_start_info, by = "read_id") 

  df <- dplyr::distinct(df)
  
  # make factor columns for categoical variables
  df %<>% dplyr::mutate(barcode = as.character(dataset)) %>% 
    dplyr::mutate(dataset = as.factor(dataset)) 
    
  df$barcode[df$barcode == '10x'] <- 10
  df$barcode[df$barcode == '15x'] <- 15
  df$barcode[df$barcode == '30x'] <- 30
  df$barcode[df$barcode == '60x'] <- 60
  df$barcode[df$barcode == '60xN'] <- 60
  df$barcode[df$barcode == '80x'] <- 80
  df$barcode[df$barcode == '100x'] <- 100
  df %<>% dplyr::mutate(barcode = as.factor(barcode)) %>% 
    dplyr::select(-transcript_alignment_start)
  
  return(df)
}
