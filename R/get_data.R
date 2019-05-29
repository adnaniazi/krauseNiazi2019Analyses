#' Download data from GitHub release using piggyback and return it as 
#' a dataframe.
#'
#' @param release_version Which version of the data to download
#' @param file_name 
#'
#' @export
#'
get_data <- function(release_version, file_name) {
  dest_dir <- file.path(here::here(), 'data')
  piggyback::pb_download(
    file = file_name,
    repo = 'adnaniazi/krauseNiazi2019Datasets',
    tag = release_version,
    dest = dest_dir
    )
  file_path <- file.path(
    dest_dir, 
    file_name
    )
  df <- read.csv(file_path, header = TRUE, stringsAsFactors = FALSE)
  return(df)
}