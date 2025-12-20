#' rename fastq files
#'
#' @description
#' This function rename files based on specific patterns. It is intended to rename
#' raw fastq files to fit with the pipeline.
#'
#' @param project_path Path to the project folder.
#' @param metadata A \code{data.frame} or \code{tibble} to change fastq names. The column
#' of the old names must be named \code{fastq}, the column for the new names \code{fastq_new}.
#' @param file_folder The folder containing r1 and r2 files.
#' @param r1_original The pattern to replace for r1.
#' @param r2_original The pattern to replace for r2.
#' @param r1_new The new pattern for r1.
#' @param r2_new The new pattern for r2.
#'
#' @details
#' The pipeline needs input fastqs end with the pattern r1.fastq. The function
#' `create_project()`try to do this for input raw files, that must be fastq.
#'
#' @export


file_rename <-  function(project_path = NULL,
                         metadata = NULL,
                         file_folder = NULL,
                         r1_original = NULL,
                         r2_original = NULL,
                         r1_new = "_r1.fastq",
                         r2_new = "_r2.fastq"){

  # set up the folder
  folder <- if (is.null(file_folder)) "1_demultiplexed" else file_folder
  folder_path <- file.path(project_path, folder)

  if (is.null(metadata)) {

    if (is.null(r1_original) || is.null(r2_original)) {
      stop("When metadata is NULL, r1_original and r2_original must be provided.")
    }

    r1_files_list <- list.files(folder_path, pattern = r1_original, full.names = TRUE)
    r2_files_list <- list.files(folder_path, pattern = r2_original, full.names = TRUE)

    r1_file_new <- gsub(r1_original, r1_new, r1_files_list)
    r2_file_new <- gsub(r2_original, r2_new, r2_files_list)

    ok1 <- file.rename(from = r1_files_list, to = r1_file_new)
    ok2 <- file.rename(from = r2_files_list, to = r2_file_new)

    if (any(!ok1) || any(!ok2)) warning("Some files could not be renamed.")
    invisible(NULL)

  } else {

    files_list <- list.files(folder_path, full.names = TRUE)
    files_list_new <- files_list

    for (i in seq_len(nrow(metadata))) {
      old_name <- metadata %>% dplyr::slice(i) %>% dplyr::pull("fastq")
      new_name <- metadata %>% dplyr::slice(i) %>% dplyr::pull("fastq_new")
      files_list_new <- gsub(old_name, new_name, files_list_new)
    }

    ok <- file.rename(from = files_list, to = files_list_new)
    if (any(!ok)) warning("Some files could not be renamed.")
    message("done")
    invisible(NULL)
  }
}
