#' Filter by quality
#'
#' @description
#' This function filter by quality using the expected error (ee) approach.
#'
#' @param project_path The path to the project folder.
#' @param vsearch_path The path to \code{vsearch} folder.
#' @param max_ee From \code{vsearch}: "...discard sequences with more than the specified number of expected errors"
#' @param qmax From \code{vsearch}: "...the maximum quality score accepted when reading FASTQ files. Default to 41."
#' @param vsearch_arguments Further arguments to be passed to \code{vsearch}.
#' @param n_reads The number of reads above which a fastq is kept for further processing.
#' See \link{names_checker} for further details.
#' @param fasta_out If \code{TRUE} return a fasta otherwise a fastq.
#'
#' @export
#'
#' @importFrom dplyr left_join
#' @importFrom writexl write_xlsx
#' @importFrom ShortRead countFastq





filter_by_quality <- function(project_path = NULL,
                              vsearch_path = NULL,
                              output_path = NULL,
                              max_ee = 1,
                              qmax = 41,
                              vsearch_arguments = NULL,
                              n_reads = 0,
                              fasta_out = TRUE){

  # get the file list to trim
  fastq_trimmed <- list.files(file.path(project_path, "5_marker_length_filtering"), full.names = TRUE) %>%
    sort()
  
  # change the names
  fastq_ee <- gsub("_merged_cutadapt_trimmed.fastq", "_ee.fastq", fastq_trimmed)

  # and the destination folder
  fastq_ee <- gsub("5_marker_length_filtering", "6_quality_filtering", fastq_ee)

  # get valid FASTQ names
  to_keep <- names_checker(project_path = project_path,
                           file_folder = "5_marker_length_filtering",
                           return_fastq_names = TRUE,
                           n_reads = n_reads) %>%
    paste(., collapse = "|")

  # select cutadapt files with grepl
  fastq_trimmed <- fastq_trimmed[grepl(to_keep,  fastq_trimmed)] %>%
    sort()

  # select trimmed files with grepl
  fastq_ee <- fastq_ee[grepl(to_keep,  fastq_ee)] %>%
    sort()

  # if fasta_out is TRUE return a fasta
  if(fasta_out){
    fastq_ee <- gsub("_ee.fastq", "_ee.fasta", fastq_ee)
  }
  
  # remove previous files
  fastq_to_remove <- list.files(file.path(output_path, "6_quality_filtering"), full.names = TRUE) %>%
    sort()
  
  if(length(fastq_to_remove > 0)){
    file.remove(fastq_to_remove)
    message("existing files removed from the folder 6_quality_filtering")
  }
  
  # create a read count data.frame after trimming
  read_count_ee <- data.frame(samples_name = gsub("_ee.fastq|_ee.fasta", "", basename(fastq_ee)),
                              ee = NA)

  # set the output path if null
  if(is.null(output_path)){
    # set the output path to project_path if the last is null
    output_path <- project_path
  }
  
  
  # change the file name according to output
  if(! is.null(output_path)){
    fastq_ee <- gsub(project_path, output_path, fastq_ee)
  }
  
  
  # create a file to store the messages from cutadapt
  file.create(file.path(output_path, "log_files", "6_quality_filtering.txt"))

  # for loop to store make trimming
  for(i in 1:length(fastq_ee)){

    if(i == 1){
      write.table(paste("###",
                        gsub("_ee.fastq", "", basename(fastq_ee[i])) ,
                        paste(rep("-", 100), collapse = "")),
                  file = file.path(output_path, "log_files", "header_primers.txt"),
                  row.names = FALSE,
                  quote = FALSE,
                  col.names = FALSE)
    } else {
      write.table(paste("\n\n\n",
                        "###",
                        gsub("_ee.fastq", "",basename(fastq_ee[i])),
                        paste(rep("-", 100), collapse = "")),
                  file = file.path(output_path, "log_files", "header_primers.txt"),
                  row.names = FALSE,
                  quote = FALSE,
                  col.names = FALSE)
    }




    file.append(file.path(output_path, "log_files", "6_quality_filtering.txt"),
                file.path(output_path, "log_files", "header_primers.txt"))



    # quality filtering
    ee_cmd <- paste("-fastq_filter ",
                    fastq_trimmed[i],
                    ifelse(fasta_out, " -fastaout ", " -fastqout "),
                    fastq_ee[i],
                    " --fastq_maxee ",
                    max_ee,
                    " --fastq_qmax ",
                    qmax,
                    vsearch_arguments,
                    sep="")

    ee_stdout <- system2(vsearch_path,
                         ee_cmd,
                         stdout = TRUE)



    writeLines(ee_stdout, file.path(output_path, "log_files", "6_quality_filtering_temp.txt"))

    file.append(file.path(project_path, "log_files", "6_quality_filtering.txt"),
                file.path(project_path, "log_files", "6_quality_filtering_temp.txt"))


    # populate the read counts file
    read_count_ee[i, 1] <- gsub("_ee.fastq|_ee.fasta", "", basename(fastq_ee[i]))
    read_count_ee[i, 2] <- ShortRead::countFastq(fastq_ee[i])[1]


    # print the progress
    message(paste(gsub("_ee.fastq", "", basename(fastq_ee[i])), ": ", i, " of ", length(fastq_ee), sep = ""))

  }


  # remove temporary files to keep the log folder clean
  file.remove(file.path(output_path, "log_files", "header_primers.txt"))
  file.remove(file.path(output_path, "log_files", "6_quality_filtering_temp.txt"))

  # add read counts to the log file
  read_count_df <- readxl::read_excel(file.path(project_path, "log_files", "0_read_count.xlsx"))

  if("ee" %in% colnames(read_count_df)){
    col_to_remove <- colnames(read_count_df)[which(colnames(read_count_df) == "ee"):ncol(read_count_df)]
    
    # remove columns with dplyr
    read_count_df %>%
      dplyr::select(-dplyr::any_of(col_to_remove)) %>%
      dplyr::left_join(., read_count_ee, by = "samples_name") %>%
      writexl::write_xlsx(., file.path(output_path, "log_files", "0_read_count.xlsx"))
  } else{
    read_count_df %>%
      dplyr::left_join(., read_count_ee, by = "samples_name") %>%
      writexl::write_xlsx(., file.path(output_path, "log_files", "0_read_count.xlsx"))
    }
}
