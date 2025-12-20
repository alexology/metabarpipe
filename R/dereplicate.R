#' Dereplicate FASTQ
#'
#' @description
#' This function dereplicates FASTQ files to FASTA files.
#'
#' @param project_path The path to the project folder.
#' @param vsearch_path The path to \code{vsearch} folder.
#' @param size_out From \code{vsearch}: "...add abundance annotations to fasta headers".
#' @param size From \code{vsearch}: "Discard sequences with a post-dereplication abundance value smaller than integer".
#' @param fastq_width From \code{vsearch}: "Fasta files produced by vsearch are wrapped (sequences are written on lines of integer nucleotides, 80 by default). Set that value to zero to eliminate the wrapping".
#' @param relabel From \code{vsearch}: "Relabel sequences using the prefix string and a ticker (1, 2, 3, etc.) to construct the new headers. Use --sizeout to conserve the abundance annotations".
#' @param vsearch_arguments Further arguments to be passed to \code{vsearch}.
#' @param n_reads The number of reads above which a fastq is kept for further processing.
#' See \link{names_checker} for further details.
#'
#' @export
#'
#' @importFrom dplyr any_of pull left_join select
#' @importFrom writexl write_xlsx
#' @importFrom Biostrings readDNAStringSet
#' @importFrom readxl read_excel


dereplicate <- function(project_path = NULL,
                        vsearch_path = NULL,
                        size_out = TRUE,
                        size = 1,
                        fastq_width = 0,
                        relabel = NULL,
                        vsearch_arguments = NULL,
                        n_reads = 0){

  # get the file list to trim
  fastq_ee <- list.files(file.path(project_path, "6_quality_filtering"), full.names = TRUE) %>%
    sort()

  # change the names
  fasta_derep <- gsub("_ee.fastq|_ee.fasta", "_derep.fasta", fastq_ee)

  # and the destination folder
  fasta_derep <- gsub("6_quality_filtering", "7_dereplication", fasta_derep)

  # get valid FASTQ names
  to_keep <- names_checker(project_path = project_path,
                           file_folder = "6_quality_filtering",
                           return_fastq_names = TRUE,
                           n_reads = n_reads) %>%
    paste(collapse = "|")



  # select cutadapt files with grepl
  fastq_ee <- fastq_ee[grepl(to_keep,  fastq_ee)] %>%
    sort()

  # select trimmed files with grepl
  fasta_derep <- fasta_derep[grepl(to_keep,  fasta_derep)] %>%
    sort()

  # create a read count data.frame after trimming
  read_count_derep <- data.frame(samples_name = gsub("_derep.fasta", "", basename(fasta_derep)),
                                 derep = NA)


  # create a file to store the messages from cutadapt
  file.create(file.path(project_path, "log_files", "7_dereplication.txt"))

  # for loop to store make trimming
  for(i in 1:length(fasta_derep)){

    if(i == 1){
      write.table(paste("###",
                        gsub("_derep.fasta", "", basename(fasta_derep[i])) ,
                        paste(rep("-", 100), collapse = "")),
                  file = file.path(project_path, "log_files", "header_primers.txt"),
                  row.names = FALSE,
                  quote = FALSE,
                  col.names = FALSE)
    } else {
      write.table(paste("\n\n\n",
                        "###",
                        gsub("_derep.fasta", "",basename(fasta_derep[i])),
                        paste(rep("-", 100), collapse = "")),
                  file = file.path(project_path, "log_files", "header_primers.txt"),
                  row.names = FALSE,
                  quote = FALSE,
                  col.names = FALSE)
    }




    file.append(file.path(project_path, "log_files", "7_dereplication.txt"),
                file.path(project_path, "log_files", "header_primers.txt"))




    cmd_der <- paste("-derep_fulllength",
                     fastq_ee[i],
                     "-output",
                     fasta_derep[i],
                     if(isTRUE(size_out)) {"-sizeout"},
                     "-minuniquesize",
                     size,
                     "-fasta_width",
                     fastq_width,
                     if(!is.null(vsearch_arguments)){vsearch_arguments},
                     if(!is.null(relabel)){paste("--relabel", relabel)},
                     sep=" ")

    derep_stdout <- system2(vsearch_path,
            cmd_der,
            stdout = TRUE)



    writeLines(derep_stdout, file.path(project_path, "log_files", "7_dereplication_temp.txt"))

    file.append(file.path(project_path, "log_files", "7_dereplication.txt"),
                file.path(project_path, "log_files", "7_dereplication_temp.txt"))


    # populate the read counts file
    read_count_derep[i, 1] <- gsub("_derep.fasta", "", basename(fasta_derep[i]))
    read_count_derep[i, 2] <- Biostrings::readDNAStringSet(fasta_derep[i]) %>%
      names() %>%
      strsplit("=", fixed = TRUE) %>%
      purrr::map_chr(2) %>%
      as.numeric() %>%
      sum()


    # print the progress
    message(paste(gsub("_derep.fasta", "", basename(fasta_derep[i])), ": ", i, " of ", length(fasta_derep), sep = ""))

  }


  # remove temporary files to keep the log folder clean
  file.remove(file.path(project_path, "log_files", "header_primers.txt"))
  file.remove(file.path(project_path, "log_files", "7_dereplication_temp.txt"))

  # add read counts to the log file
  read_count_df <- readxl::read_excel(file.path(project_path, "log_files", "0_read_count.xlsx"))

  if("derep" %in% colnames(read_count_df)){
    col_to_remove <- colnames(read_count_df)[which(colnames(read_count_df) == "derep"):ncol(read_count_df)]
    
    # remove columns with dplyr
    read_count_df %>%
      dplyr::select(-dplyr::any_of(col_to_remove)) %>%
      dplyr::left_join(read_count_derep, by = "samples_name") %>%
      writexl::write_xlsx(file.path(project_path, "log_files", "0_read_count.xlsx"))
  } else{
    read_count_df %>%
      dplyr::left_join(read_count_derep, by = "samples_name") %>%
      writexl::write_xlsx(file.path(project_path, "log_files", "0_read_count.xlsx"))
  }

}
