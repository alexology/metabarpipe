#' Display starting checks
#'
#' @description
#' This function display starting checks in the viewer pane.
#'
#' @param project_path Path to the project folder.
#'
#' @export
#'
#' @importFrom readxl read_excel
#' @importFrom knitr kable
#' @importFrom kableExtra kable_styling

display_starting_checks <- function(project_path = NULL){

  # print the result to viewer
  readxl::read_excel(file.path(project_path, "log_files", "1_initial_checks.xlsx")) %>%
    knitr::kable() %>%
    kableExtra::kable_styling("striped")
}
