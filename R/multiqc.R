#' Run multiqc
#'
#' @description
#' This function run the python progam multiqc.
#'
#' @details
#' The function `multiqc` run the python program multiqc. Python and multiqc must
#' be already installed. See \url{https://multiqc.info/} for more details.
#'
#' @param project_path Path to the project folder.
#'
#' @returns
#' The function `multiqc`return an html file that is open in the default browser
#'
#' @export
#'
#' @importFrom utils read.table

multiqc <- function(project_path = NULL){

  # read the project path
  project_path <- utils::read.table(file.path(project_path, "project_path.txt")) %>%
    as.character()

  # get the current working directory
  current_wd <- getwd()

  # set the working directory to avoid problems with multiqc
  setwd(file.path(project_path, "reports"))

  # run multiqc
  system("py -m multiqc . ")

  # open the html from R
  browseURL(paste0('file://', file.path(project_path, "reports", "multiqc_report.html")))

  # restore the working directory
  setwd(current_wd)

}
