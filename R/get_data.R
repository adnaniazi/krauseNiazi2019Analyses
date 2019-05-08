#' Download data from GitHub release using piggyback.
#'
#' @param release_version Which version of the data to download
#'
#' @export
#'
get_data <- function(release_version, file_name) {
  dest_dir <- system.file('extdata', package = 'krauseNiazi2019Analyses')
  piggyback::pb_download(
    file = file_name,
    repo = 'adnaniazi/krauseNiazi2019Datasets',
    tag = release_version,
    dest = dest_dir
    )
  file_path <- system.file(
    'extdata', 
    file_name,
    package = 'krauseNiazi2019Analyses')
  df <- read.csv(file_path, header = TRUE, stringsAsFactors = FALSE)
}