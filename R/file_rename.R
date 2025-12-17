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
#'
#' @importFrom dplyr select slice pull

file_rename <-  function(project_path = NULL,
                         metadata = NULL,
                         file_folder = NULL,
                         r1_original = NULL,
                         r2_original = NULL,
                         r1_new = "_r1.fastq",
                         r2_new = "_r2.fastq"){

  if(is.null(metadata)){
    # get r1 file list
    r1_files_list <- list.files(file.path(project_path, "1_demultiplexed"),
                                pattern = r1_original,
                                full.names = TRUE)

    # get r2 file list
    r2_files_list <- list.files(file.path(project_path, "1_demultiplexed"),
                                pattern = r2_original,
                                full.names = TRUE)

    # change the pattern of r1
    r1_file_new <- gsub(r1_original, r1_new, r1_files_list)

    # change the pattern of r2
    r2_file_new <- gsub(r2_original, r2_new, r2_files_list)

    # rename r1
    file.rename(from = r1_files_list, to = r1_file_new)

    # rename r2
    file.rename(from = r2_files_list, to = r2_file_new)


  } else {
    # get file list
    files_list <- list.files(file.path(project_path, file_folder),
                             full.names = TRUE)

    # copy the old file list to generate the new one
    files_list_new <- files_list

    for(i in 1:nrow(metadata)){

      # get the old name
      old_name <-  metadata %>%
        dplyr::select(fastq) %>%
        dplyr::slice(i) %>%
        dplyr::pull()


      # get new name
      new_name <-  metadata %>%
        dplyr::select(fastq_new) %>%
        dplyr::slice(i) %>%
        dplyr::pull()


      # change the names
      files_list_new <- gsub(old_name, new_name, files_list_new)

    }

    file.rename(from = files_list, to = files_list_new)
    message("done")
  }




}
