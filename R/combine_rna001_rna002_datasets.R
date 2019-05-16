#' Title
#'
#' @param rna001_dataset 
#' @param rna002_dataset 
#'
#' @return
#' @export
#'
#' @examples
combine_rna001_rna002_datasets <- function(rna001_dataset,
                                           rna002_dataset) {
  rna001_dataset %<>% 
    # add a kit column: both replicates 1 and 2 in this dataset were
    # sequenced using RNA001 kit
    dplyr::mutate(kit = 'SQK-RNA001')
  
  rna002_dataset %<>% 
    # add a kit column
    # It was sequenced using RNA002 kit
    dplyr::mutate(kit = 'SQK-RNA002') %>% 
    # add a replicate column
    # technically this is the third replicate
    dplyr::mutate(replicate = 3) 
    
  # bind both these datasets
  df <- rbind(rna001_dataset, rna002_dataset)
  
  # make factor columns for categoical variables
  df %<>% dplyr::mutate(replicate = as.factor(replicate), 
                        barcode = as.factor(barcode),
                        kit = as.factor(kit))
  
  return(df)
}