#' Check fastq names from r1 and r2
#'
#' @description
#' This function checks if r1 and r2 files match.
#'
#' @param project_path Path to the project folder.
#' @param file_folder The folder containing r1 and r2 files.
#' @param r1_pattern Pattern to remove from r1 to make the comparison with r2.
#' Default to `_r1.fastq.`
#' @param r2_pattern Pattern to remove from r2 to make the comparison with r1.
#' Default to `_r2.fastq.`
#' @param return_fastq_names if \code{TRUE}, returns a list of fastq files present
#' in both r1 and r2.
#' @param n_reads The number of reads above which a fastq is kept for further processing.
#'
#' @details
#' The function `r1_r2_checker` remove patterns from r1 and r2 to check if the
#' file names from r1 and r2 are the same.
#'
#' @export
#'
#' @importFrom dplyr filter pull rename
#' @importFrom readxl read_excel
#'
#' @returns
#' The function `r1_r2_checker` returns TRUE if all the r1 files are matched by
#' files in r2, otherwise FALSE. If \code{return_fastq_names} is set to \code{TRUE},
#' the function returns a list of r1 and r2 fastq.


names_checker <- function(project_path = NULL,
                          file_folder = NULL,
                          r1_pattern = "_r1.fastq",
                          r2_pattern = "_r2.fastq",
                          return_fastq_names = FALSE,
                          n_reads = 0){

  if(file_folder %in% c("1_demultiplexed", "2_trim_by_length")){
    # change the pattern of r1
    r1_files <- list.files(file.path(project_path, file_folder),
                           pattern = r1_pattern,
                           full.names = FALSE) %>%
      sort()

    # change the pattern of r2
    r2_files <- list.files(file.path(project_path, file_folder),
                           pattern = r2_pattern,
                           full.names = FALSE) %>%
      sort()

    # remove the pattern from r1 to make the names comparable
    r1_files_list <- gsub(r1_pattern, "", r1_files)

    # remove the pattern from r2 to make the names comparable
    r2_files_list <- gsub(r2_pattern, "", r2_files)

    if(!isTRUE(return_fastq_names)){
      # check if files are the same
      if(identical(r1_files_list, r2_files_list)){
        return(TRUE)
      } else {
        return(FALSE)
      }
    }

    if(isTRUE(return_fastq_names)){
      read_count_df <- readxl::read_excel(file.path(project_path, "log_files", "0_read_count.xlsx"))
      reads_intersect <- intersect(r1_files_list, r2_files_list)


      to_subset <- switch(file_folder,
                          "1_demultiplexed" = c("samples_name", "r1_raw", "r2_raw"),
                          "2_trim_by_length" = c("samples_name", "r1_trimmed", "r2_trimmed"))


      to_keep <- read_count_df[, to_subset] %>%
          na.omit() %>%
          dplyr::rename(r1 = 2, r2 = 3) %>%
          dplyr::filter(samples_name %in% reads_intersect) %>%
          dplyr::filter((r1 > n_reads) & (r2 > n_reads)) %>%
          pull(samples_name)


      if(length(to_keep) == 0){
        stop("no file selected")
      } else {

        # returns the name to keep
        return(to_keep)

      }




    }
  }



  if(file_folder %in% c("3_paired_end_merge",
                        "4_cutadapt_trimmed",
                        "5_marker_length_filtering",
                        "6_quality_filtering")){

    fastq_files <- list.files(file.path(project_path, file_folder),
                              full.names = FALSE)

    read_count_df <- readxl::read_excel(file.path(project_path, "log_files", "0_read_count.xlsx"))

    to_subset <- switch(file_folder,
                        "3_paired_end_merge" = c("samples_name", "merged"),
                        "4_cutadapt_trimmed" = c("samples_name", "cutadapt"),
                        "5_marker_length_filtering" = c("samples_name", "marker_trim"),
                        "6_quality_filtering" = c("samples_name", "ee"))


    to_keep <- read_count_df[, to_subset] %>%
      na.omit() %>%
      dplyr::rename(r1 = 2) %>%
      dplyr::filter(r1 > n_reads) %>%
      dplyr::pull(samples_name)


    if(length(to_keep) == 0){
      stop("no file selected")
    } else {

      # returns the names to keep

      return(to_keep)


    }

  }

}
