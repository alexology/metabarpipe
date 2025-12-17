#' count fastq reads
#'
#' @description
#' This function open the log file about read counts stored in the log_files folder.
#'
#' @param project_path Path to the project folder.
#'
#' @export
#'
#' @importFrom readxl read_excel
#' @importFrom knitr kable
#' @importFrom kableExtra kable_styling

count_reads <-  function(project_path = NULL){

  # print the result to viewer
  readxl::read_excel(file.path(project_path, "log_files", "0_read_count.xlsx")) %>%
    knitr::kable() %>%
    kableExtra::kable_styling("striped")
}
