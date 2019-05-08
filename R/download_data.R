#' Download data from GitHub release using piggyback.
#'
#' @param release_version Which version of the data to download
#'
#' @export
#'
download_data <- function(release_version) {
  dest_dir <- system.file('extdata', package = 'krauseNiazi2019Analyses')
  piggyback::pb_download(repo = 'adnaniazi/krauseNiazi2019Datasets',
                         tag = release_version,
                         dest = dest_dir)
  
}