#' Title
#'
#' @param decoded_barcodes 
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
consolidate_krause_etal_rna001_data <- function(decoded_barcodes,
                                                tailfindr_estimates,
                                                nanopolish_estimates,
                                                transcript_start_info) {
  
  decoded_barcodes <- decoded_barcodes %>%
    # keep only valid reads in which GFP sequence is there
    dplyr::filter(read_type == 'polyA') %>% 
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
    dplyr::select(read_id, barcode, replicate, file_path)
    
  tailfindr_estimates <- tailfindr_estimates %>% 
    # filter out reads for which tailfindr's tail length
    # estimates are not available
    dplyr::filter(!is.na(tail_length)) %>% 
    # remove irrelevant columns
    dplyr::select(-polya_fastq, -file_path, -replicate) %>% 
    # rename variables:
    dplyr::rename(tail_start_tf = tail_start,
                  tail_end_tf = tail_end,
                  samples_per_nt_tf = samples_per_nt,
                  tail_length_tf = tail_length)
    
  nanopolish_estimates <- nanopolish_estimates %>% 
    # # for read which have nanopolish's prediction as -1, set
    # # these predictions to NA
    # dplyr::mutate(polya_length = ifelse(polya_length == -1, NA,  polya_length),
    #               read_rate = ifelse(read_rate == -1, NA,  read_rate),
    #               polya_start = ifelse(polya_start == -1, NA,  polya_start),
    #               transcript_start = ifelse(transcript_start == -1, NA,  transcript_start)) %>%
    # take only PASS reads
    dplyr::filter(qc_tag == 'PASS') %>% 
    # reomve irrelevant columns
    dplyr::select(-contig, -position, -leader_start, 
                  -adapter_start, -replicate, -file_path) %>% 
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
  
  # join all these four datasets
  # we don't need tailfindr's or nanopolish's predication for reads
  # that we can't find any barcodes for. So merge by barcode infomation first
  df <- dplyr::left_join(decoded_barcodes, tailfindr_estimates, by = "read_id") %>% 
    dplyr::left_join(nanopolish_estimates, by = "read_id") %>% 
    dplyr::left_join(transcript_start_info, by = "read_id") 
  
  df <- dplyr::distinct(df)
  
  return(df)
}
