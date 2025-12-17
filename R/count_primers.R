#' count primer occurrence
#'
#' @description
#' This function prints the number of occurrence in fastq files to the viewer pane.
#'
#' @param project_path Path to the project folder.
#' @param orientation The orientation. Default to r1.
#'
#' @details
#' The information about the number of occurrence is stored in the log_files folder.
#'
#'
#' @export
#'
#' @importFrom utils read.table
#' @importFrom readxl read_excel
#' @importFrom knitr kable
#' @importFrom kableExtra kable_styling

count_primers <- function(project_path = NULL,
                          orientation = "r1"){

  if(identical(orientation, "r1")){
    # read the project path
    project_path <- read.table(file.path(project_path, "project_path.txt")) %>%
      as.character()

    # print the result to viewer
    readxl::read_excel(file.path(project_path, "log_files", "2_primer_count_original_fastq_R1.xlsx")) %>%
      knitr::kable() %>%
      kableExtra::kable_styling("striped")
  } else {
    # read the project path
    project_path <- read.table(file.path(project_path, "project_path.txt")) %>%
      as.character()

    # print the result to viewer
    readxl::read_excel(file.path(project_path, "log_files", "2_primer_count_original_fastq_R2.xlsx")) %>%
      knitr::kable() %>%
      kableExtra::kable_styling("striped")
  }


}
