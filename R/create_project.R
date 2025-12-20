#' create a project
#'
#' @description
#' This function creates a project to run the metabarcoding pipeline.
#'
#' @param project_path The path to the project folder. The folder must already exists.
#'
#' @details
#' The function `create_project()` creates 12 folders in the project folder.
#' A txt file storing the project path is also created in the log_files folder.
#' The path should be absolute rather than relative.
#'
#' @importFrom utils write.table
#'
#' @export
#'
#' @examples
#' \dontrun{
#' my_folder_path <- "D:/my_project"
#' create_project(my_folder_path)
#' }


create_project <- function(project_path){
  folders <- c("log_files",
               "reports",
               "0_raw_files",
               "1_demultiplexed",
               "2_trim_by_length",
               "3_paired_end_merge",
               "4_cutadapt_trimmed",
               "5_marker_length_filtering",
               "6_quality_filtering",
               "7_dereplication",
               "8_denoising",
               "9_chimera",
               "10_asv_table",
               "11_translation",
               "12_results",
               "13_taxonomic_assignment")

  # create folders
  lapply(folders, function(x) dir.create(file.path(project_path, x), showWarnings = FALSE))

  # write file that store the log path
  utils::write.table(project_path,
                     file = file.path(project_path, "project_path.txt"),
                     col.names = FALSE,
                     row.names = FALSE)

  print("Done")
}
