#' Title
#'
#' @param event_data 
#' @param start 
#' @param end 
#'
#' @return
#' @export
#'
#' @examples
get_total_moves_between_borders <- function(event_data, start, end) {
    row_index_start <- which.min(abs(event_data$start - start))
    row_index_end <- which.min(abs(event_data$start - end))
    total_moves <- sum(event_data$move[row_index_start:row_index_end])
    return(total_moves)
}