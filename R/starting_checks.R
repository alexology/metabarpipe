#' Further quality checks
#'
#' @description
#' This function compare the number of reads and the reads name between
#' r1 and r2
#'
#' @param project_path Path to the project folder.
#' @param r1_pattern The pattern of r1. Default to \code{_r1.fastq}.
#' @param r2_pattern The pattern of r1. Default to \code{_r2.fastq}.
#' @param n The number of reads to compare the read names of r1 and r2.
#'
#' @details
#' The function \code{starting checks} compare the read names of r1 and r2 on a
#' subset of the reads. Increase \code{n} if you need more accurate results.
#'
#' @export
#'
#' @importFrom writexl write_xlsx
#' @importFrom ShortRead countFastq readFastq
#' @importFrom dplyr select


starting_checks <- function(project_path = NULL,
                            r1_pattern = "_r1.fastq",
                            r2_pattern = "_r2.fastq",
                            n = 10e3){

  # get r1 file list
  r1_files_list <- list.files(path = file.path(project_path, "1_demultiplexed"),
                              pattern = r1_pattern,
                              full.names = TRUE) %>%
    sort()

  # get r2 file list
  r2_files_list <- list.files(path = file.path(project_path, "1_demultiplexed"),
                              pattern = r2_pattern,
                              full.names = TRUE) %>%
    sort()

  # get the samples name
  sample_names_r1 <- gsub("_r1.fastq", "", basename(r1_files_list))


  # create a data.frame to store the results
  # length_check will check if the fastq length of r1 and r2 are the same, TRUE if
  # they are the same
  # id_check will do the same for the ID, TRUE if they are the same
  initial_checks <- data.frame(samples_name = sample_names_r1,
                               length_check = NA,
                               id_check = NA)

  # create a read count data.frame
  read_count <- data.frame(samples_name = sample_names_r1,
                           r1_raw = NA,
                           r2_raw = NA)

  # iterate over all the files
  for(i in 1:length(r1_files_list)){

    # get fastq length for r1 and r2
    temp_length_r1 <- ShortRead::countFastq(r1_files_list[i])
    temp_length_r2 <- ShortRead::countFastq(r2_files_list[i])

    # check if r1 and r2 have the same length
    initial_checks[i, "length_check"] <- temp_length_r1[1] == temp_length_r2[1]

    # store the read counts
    read_count[i, "r1_raw"] <- temp_length_r1[1]
    read_count[i, "r2_raw"] <- temp_length_r2[1]

    # take subset of r1 and r2 to check for IDs
    subset_r1 <- ShortRead::readFastq(r1_files_list[i], n = n)
    subset_r2 <- ShortRead::readFastq(r2_files_list[i], n = n)

    # get the IDs without the last part for r1
    id_r1_temp <- subset_r1@id %>%
      as.character() %>%
      strsplit(., " ") %>%
      do.call(rbind, .) %>%
      as.data.frame() %>%
      dplyr::select(1)

    # get the IDs without the last part for r2
    id_r2_temp <- subset_r2@id  %>%
      as.character() %>%
      strsplit(., " ") %>%
      do.call(rbind, .) %>%
      as.data.frame() %>%
      dplyr::select(1)

    # check if IDs are the same
    initial_checks[i, "id_check"] <- all(id_r1_temp == id_r2_temp)

    # print the progress
    print(paste(sample_names_r1[i], ": ", i, " of ", length(sample_names_r1), sep = ""))

  }

  # test if checks failed
  if(all(initial_checks[, 2:3])){"check passed"} else {warning("check failed")}

  # write check results to disk
  writexl::write_xlsx(initial_checks, file.path(project_path, "log_files", "1_initial_checks.xlsx"))

  # write the read count
  writexl::write_xlsx(read_count, file.path(project_path, "log_files", "0_read_count.xlsx"))
}
