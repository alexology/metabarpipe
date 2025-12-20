#' Trim by length
#'
#' @description
#' This function trim the left or the right side of the reads in the raw fastq.
#'
#' @param project_path The path to the project folder.
#' @param vsearch_path The path to \code{vsearch} folder.
#' @param usearch_path The path to \code{usearch} folder.
#' @param n_reads The number of reads above which a fastq is kept for further processing.
#' @param trim_r1 Trim length of r1.
#' @param trim_r2 Trim length of r2.
#' @param left_r1 Number of bases to remove on the left r1.
#' @param right_r1 Number of bases to remove on the right r1.
#' @param left_r2 Number of bases to remove on the left r2.
#' @param right_r2 Number of bases to remove on the right r2.
#' @param min_length_r1 Minum length of the reads to keep for r1.
#' @param min_length_r2 Minum length of the reads to keep for r2.
#' @param vsearch_arguments_r1 Further arguments to be passed to \code{vsearch} for r1.
#' @param vsearch_arguments_r2 Further arguments to be passed to \code{vsearch} for r2.
#'
#' @details
#' The function \code{trim_fastq_length} relies on \code{vsearch}. It is intended
#' for trimming raw fastq files. Sometimes the reads from one orientation can be
#' trimmed more than the other, resulting in a mismatch between reads ID of the
#' two orientations. When a mismatch is detected, the function \code{trim_fastq_length}
#' perform a simultaneus trimming between the two orientation. With the simultaneous option
#' do not allow different trimming between r1 and r2.  This behavior will likely
#' change in the future.
#'
#' @export
#'
#' @importFrom ShortRead countFastq
#' @importFrom Biostrings readDNAStringSet width
#' @importFrom dplyr select left_join
#' @importFrom readxl read_excel
#' @importFrom writexl write_xlsx




trim_fastq_length <- function(project_path = NULL,
                              vsearch_path = NULL,
                              usearch_path = NULL,
                              n_reads = 0,
                              trim_r1 = NULL,
                              trim_r2 = NULL,
                              left_r1 = 0,
                              right_r1 = 0,
                              left_r2 = 0,
                              right_r2 = 0,
                              min_length_r1 = NULL,
                              min_length_r2 = NULL,
                              vsearch_arguments_r1 = NULL,
                              vsearch_arguments_r2 = NULL){


  # get the path to r1 files
  r1_path <- list.files(file.path(project_path, "1_demultiplexed"),
                        full.names = TRUE,
                        pattern = "_r1.fastq") %>%
    sort()

  # get the path to r2 files
  r2_path <- list.files(file.path(project_path, "1_demultiplexed"),
                        full.names = TRUE,
                        pattern = "_r2.fastq") %>%
    sort()

  # select the files based on names_checker
  to_keep <- names_checker(project_path = project_path,
                           file_folder = "1_demultiplexed",
                           return_fastq_names = TRUE,
                           n_reads = n_reads) %>%
    paste(collapse = "|")


  # select r1 files with grepl
  r1_path <- r1_path[grepl(to_keep,  r1_path)] %>%
    sort()

  # select r2 files with grepl
  r2_path <- r2_path[grepl(to_keep,  r2_path)] %>%
    sort()

  # path for r1 trimmed data
  r1_new_path <- gsub("_r1.fastq", "_r1_trimmed.fastq", r1_path)
  r1_new_path <- gsub("1_demultiplexed", "2_trim_by_length", r1_new_path)

  # path for r2 trimmed data
  r2_new_path <- gsub("_r2.fastq", "_r2_trimmed.fastq", r2_path)
  r2_new_path <- gsub("1_demultiplexed", "2_trim_by_length", r2_new_path)

  # set a data.frame to store count results
  read_count_trimmed <- data.frame(samples_name = NA,
                                   r1_trimmed = NA,
                                   r2_trimmed = NA)


  for(i in 1:length(r1_path)){
    # set the vsearch commands for r1
    cmd_r1 <- paste("--fastx_filter ",
                    r1_path[i],
                    " --fastq_stripleft ", left_r1,
                    " --fastq_stripright ", right_r1,
                    if(!is.null(trim_r1)){paste(" --fastq_trunclen ", trim_r1, sep = "")},
                    # if(is.null(min_length_r1)){paste(" --fastq_minlen ", min_length_r1, sep = "")},
                    if(!is.null(vsearch_arguments_r1)){paste(" ", vsearch_arguments_r1, " ", sep = "")},
                    " --fastqout ",
                    r1_new_path[i],
                    sep="")

    # set the vsearch commands for r2
    cmd_r2 <- paste("--fastx_filter ",
                    r2_path[i],
                    " --fastq_stripleft ", left_r2,
                    " --fastq_stripright ", right_r2,
                    if(!is.null(trim_r2)){paste(" --fastq_trunclen ", trim_r2, sep = "")},
                    # if(is.null(min_length_r2)){paste(" --fastq_minlen ", min_length_r2, sep = "")},
                    if(!is.null(vsearch_arguments_r2)){paste(" ", vsearch_arguments_r2, " ", sep = "")},
                    " --fastqout ",
                    r2_new_path[i],
                    sep="")


    # run vsearch for r1
    s2_r1 <- system2(vsearch_path,
                     cmd_r1)

    # run vsearch for r2
    s2_r2 <- system2(vsearch_path,
                     cmd_r2)

    # take subset of r1 and r2 to check for IDs
    subset_r1 <- Biostrings::readDNAStringSet(r1_new_path[i],
                                              nrec = 10e2,
                                              format = "fastq")
    subset_r2 <- Biostrings::readDNAStringSet(r2_new_path[i],
                                              nrec = 10e2,
                                              format = "fastq")

    # get the IDs without the last part for r1
    id_r1_temp <- names(subset_r1) %>%
      as.character() %>%
      strsplit(" ")
    
    # to avoid RCMD checks
    id_r1_temp <- do.call(rbind, id_r1_temp) %>%
      as.data.frame() %>%
      dplyr::select(1)
    

    # get the IDs without the last part for r2
    id_r2_temp <- names(subset_r2)  %>%
      as.character() %>%
      strsplit(" ")
    
    # to avoid RCMD checks
    id_r2_temp <-  do.call(rbind, id_r2_temp) %>%
      as.data.frame() %>%
      dplyr::select(1)

    r1_read_count <- ShortRead::countFastq(r1_new_path[i])[1]
    r2_read_count <- ShortRead::countFastq(r2_new_path[i])[1]

    if( (! all(id_r1_temp == id_r2_temp)) | (r1_read_count != r2_read_count)){
      message("some mess with IDs")
      closeAllConnections()

      file_fwd <- r1_new_path[i]
      file_rvr <- r2_new_path[i]

      file_fwd <- gsub("_r1_trimmed.fastq", "_r1_original.fastq", file_fwd)
      file_rvr <- gsub("_r2_trimmed.fastq", "_r2_original.fastq", file_rvr)

      file.rename(from = r1_new_path[i], to = file_fwd)
      file.rename(from = r2_new_path[i], to = file_rvr)

      cmd_usearch <- paste("-fastx_syncpairs",
                           file_fwd,
                           "-reverse",
                           file_rvr,
                           "-output",
                           r1_new_path[i],
                           "-output2",
                           r2_new_path[i])

      system2(usearch_path,
              cmd_usearch)

      closeAllConnections()
      file.remove(file_fwd)
      file.remove(file_rvr)



      # take subset of r1 and r2 to check for IDs
      subset_r1 <- Biostrings::readDNAStringSet(r1_new_path[i],
                                                nrec = 10e2,
                                                format = "fastq")
      subset_r2 <- Biostrings::readDNAStringSet(r2_new_path[i],
                                                nrec = 10e2,
                                                format = "fastq")

      # get the IDs without the last part for r1
      id_r1_temp <- names(subset_r1) %>%
        as.character() %>%
        strsplit(" ")
      
      # to avoid RCMD checks
      id_r1_temp <- do.call(rbind, id_r1_temp) %>%
        as.data.frame() %>%
        dplyr::select(1) %>%
        dplyr::pull(1)

      # get the IDs without the last part for r2
      id_r2_temp <- names(subset_r2)  %>%
        as.character() %>%
        strsplit(" ")
      
      # to avoid RCMD checks
      id_r2_temp <- do.call(rbind, id_r2_temp) %>%
        as.data.frame() %>%
        dplyr::select(1) %>%
        dplyr::pull(1)

      r1_read_count <- ShortRead::countFastq(r1_new_path[i])[1]
      r2_read_count <- ShortRead::countFastq(r2_new_path[i])[1]


      if(all(id_r1_temp == id_r2_temp)){
        message("mess fixed")
      } else{
        message("still mess")
      }
    }


    # check if trimming length is right
    if(max(Biostrings::width(subset_r1)) > trim_r1){
      message("Something wrong with the truncation. Please check the code.")
    }

    if(max(Biostrings::width(subset_r2)) > trim_r2){
      message("Something wrong with the truncation. Please check the code.")
    }

    # remove objects to save space
    rm(subset_r1)
    rm(subset_r2)

    # populate the read counts file
    read_count_trimmed[i, 1] <- gsub("_r1_trimmed.fastq", "", basename(r1_new_path[i]))
    read_count_trimmed[i, 2] <- r1_read_count
    read_count_trimmed[i, 3] <- r2_read_count

    message(paste(gsub("_r1_trimmed.fastq", "", basename(r1_new_path[i])), ": ", i, " of ", length(r1_new_path), sep = ""))


  }


  read_count_df <- readxl::read_excel(file.path(project_path, "log_files", "0_read_count.xlsx"))

  if(any(c("r1_trimmed", "r2_trimmed") %in% colnames(read_count_df))){
    read_count_df[, "r1_trimmed"] <- read_count_trimmed[,"r1_trimmed"]
    read_count_df[, "r2_trimmed"] <- read_count_trimmed[,"r2_trimmed"]
    writexl::write_xlsx(read_count_df, file.path(project_path, "log_files", "0_read_count.xlsx"))
  } else{
    read_count_df  %>%
      dplyr::left_join(read_count_trimmed, by = "samples_name") %>%
      writexl::write_xlsx(file.path(project_path, "log_files", "0_read_count.xlsx"))
  }

  closeAllConnections()
}
