#' Title
#'
#' @param read_path 
#' @param plot_debug 
#' @param basecalled_with 
#' @param multifast5 
#' @param model 
#' @param plotting_library 
#'
#' @return
#' @export
#'
#' @examples
get_moves_per_read <- function(file_path,
                               basecalled_with,
                               multifast5,
                               model,
                               tail_start,
                               tail_end) {
    
    read_data <- extract_read_data(file_path = file_path,
                                   read_id_fast5_file = NA,
                                   plot_debug = F,
                                   basecalled_with = basecalled_with,
                                   basecall_group = 'Basecall_1D_000',
                                   multifast5 = multifast5,
                                   model = model,
                                   plotting_library = 'ggplot2') 
    
    event_data <- read_data$event_data
    moves_in_tail <- get_total_moves_between_borders(event_data, tail_start, tail_end)
    read_id <- read_data$read_id
    
    list(read_id = read_id,
         moves_in_tail = moves_in_tail)
}