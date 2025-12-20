#' Trim to marker length
#'
#' @description
#' Thi function trim reads to marker length.
#'
#' @param project_path The path to the project folder.
#' @param cutadapt_path The path to \code{cutadapt} folder.
#' @param m Minimum length.
#' @param M Maximum length.
#' @param n_cores The number of cores.
#' @param n_reads The number of reads above which a fastq is kept for further processing.
#' See \link{names_checker} for further details.
#'
#' @export
#'
#' @importFrom ShortRead countFastq
#' @importFrom dplyr left_join
#' @importFrom writexl write_xlsx


trim_to_marker_length <- function(project_path = NULL,
                                  cutadapt_path = NULL,
                                  m = NULL,
                                  M = NULL,
                                  n_cores = 1,
                                  n_reads = 0){

  # get the file list to trim
  fastq_cutadapt <- sort(list.files(file.path(project_path, "4_cutadapt_trimmed"), full.names = TRUE))

  # change the names
  fastq_trimmed <- gsub("_merged_cutadapt.fastq", "_merged_cutadapt_trimmed.fastq", fastq_cutadapt)

  # and the destination folder
  fastq_trimmed <- gsub("4_cutadapt_trimmed", "5_marker_length_filtering", fastq_trimmed)


  to_keep <- names_checker(project_path = project_path,
                           file_folder = "4_cutadapt_trimmed",
                           return_fastq_names = TRUE,
                           n_reads = n_reads) %>%
    paste(collapse = "|")


  # select cutadapt files with grepl
  fastq_cutadapt <- fastq_cutadapt[grepl(to_keep,  fastq_cutadapt)] %>%
    sort()

  # select trimmed files with grepl
  fastq_trimmed <- fastq_trimmed[grepl(to_keep,  fastq_trimmed)] %>%
    sort()

  # create a read count data.frame after trimming
  read_marker_trim <- data.frame(samples_name = gsub("_merged_cutadapt.fastq", "", basename(fastq_trimmed)),
                                 marker_trim = NA)


  # create a file to store the messages from cutadapt
  file.create(file.path(project_path, "log_files", "5_trim_to_marker_length.txt"))

  # for loop to store make trimming
  for(i in 1:length(fastq_trimmed)){

    if(i == 1){
      write.table(paste("###",
                        gsub("_merged_cutadapt_trimmed.fastq", "", basename(fastq_trimmed[i])) ,
                        paste(rep("-", 100), collapse = "")),
                  file = file.path(project_path, "log_files", "header_primers.txt"),
                  row.names = FALSE,
                  quote = FALSE,
                  col.names = FALSE)
    } else {
      write.table(paste("\n\n\n",
                        "###",
                        gsub("_merged_cutadapt_trimmed.fastq", "",basename(fastq_trimmed[i])),
                        paste(rep("-", 100), collapse = "")),
                  file = file.path(project_path, "log_files", "header_primers.txt"),
                  row.names = FALSE,
                  quote = FALSE,
                  col.names = FALSE)
    }

    file.append(file.path(project_path, "log_files", "5_trim_to_marker_length.txt"),
                file.path(project_path, "log_files", "header_primers.txt"))

    # compose cutadapt command
    trim_cmd <- paste(fastq_cutadapt[i],
                      "-o ", fastq_trimmed[i],
                      "-m ", m,
                      "-M ", M,
                      "-j", n_cores)

    # run the command
    trim_stdout <- system2(cutadapt_path,
            args = trim_cmd,
            stdout = TRUE)


    trim_stdout <- trim_stdout[suppressWarnings(which((grepl("=== Summary ===", trim_stdout))):length(trim_stdout))]

    writeLines(trim_stdout, file.path(project_path, "log_files", "5_trim_to_marker_length_temp.txt"))

    file.append(file.path(project_path, "log_files", "5_trim_to_marker_length.txt"),
                file.path(project_path, "log_files", "5_trim_to_marker_length_temp.txt"))


    # populate the read counts file
    read_marker_trim[i, 1] <- gsub("_merged_cutadapt_trimmed.fastq", "", basename(fastq_trimmed[i]))
    read_marker_trim[i, 2] <- ShortRead::countFastq(fastq_trimmed[i])[1]



    # print the progress
    message(paste(gsub("_merged_cutadapt_trimmed.fastq", "", basename(fastq_trimmed[i])), ": ", i, " of ", length(fastq_trimmed), sep = ""))

  }


  # remove temporary files to keep the log folder clean
  file.remove(file.path(project_path, "log_files", "header_primers.txt"))
  file.remove(file.path(project_path, "log_files", "5_trim_to_marker_length_temp.txt"))


  # add read counts to the log file
  read_count_df <- readxl::read_excel(file.path(project_path, "log_files", "0_read_count.xlsx"))

    if("marker_trim" %in% colnames(read_count_df)){
    col_to_remove <- colnames(read_count_df)[which(colnames(read_count_df) == "marker_trim"):ncol(read_count_df)]
    
    # remove columns with dplyr
    read_count_df %>%
      dplyr::select(-dplyr::any_of(col_to_remove)) %>%
      dplyr::left_join(read_marker_trim, by = "samples_name") %>%
      writexl::write_xlsx(file.path(project_path, "log_files", "0_read_count.xlsx"))
  } else{
    read_count_df %>%
      dplyr::left_join(read_marker_trim, by = "samples_name") %>%
      writexl::write_xlsx(file.path(project_path, "log_files", "0_read_count.xlsx"))
  }

  closeAllConnections()
}
