#' Denoise
#'
#' @description
#' This function performs denoising.
#'
#' @param project_path The path to the project folder.
#' @param vsearch_path The path to \code{vsearch} folder.
#' @param size_in From \code{vsearch}: "In de novo mode, abundance annotations present in sequence headers are taken into account by default".
#' @param size_out From \code{vsearch}: "...add abundance annotations to fasta headers".
#' @param minsize From \code{vsearch}: "Specify the minimum abundance of sequences for denoising using --cluster_unoise. The default is 8".
#' @param unoise_alpha From \code{vsearch}: "Specify the alpha parameter to the --cluster_unoise command. The default i 2".
#' @param fasta_width From \code{vsearch}: "Fasta files produced by vsearch are wrapped (sequences are written on lines of integer nucleotides, 80 by default). Set that value to zero to eliminate the wrapping".
#' @param vsearch_arguments Further arguments to be passed to \code{vsearch}.
#' @param n_reads The number of reads above which a fastq is kept for further processing.
#' See \link{names_checker} for further details.
#'
#' @export
#'
#' @importFrom dplyr pull left_join
#' @importFrom writexl write_xlsx
#' @importFrom Biostrings readDNAStringSet

denoise <-  function(project_path = NULL,
                     vsearch_path = NULL,
                     size_in = TRUE,
                     size_out = TRUE,
                     minsize = 8,
                     unoise_alpha = 2,
                     fasta_width = 0,
                     vsearch_arguments = NULL,
                     n_reads = 0){

  # get the file list to trim
  fasta_derep <- list.files(file.path(project_path, "7_dereplication"), full.names = TRUE) %>%
    sort()

  # change the names
  fasta_dnoise <- gsub("_derep.fasta", "_dnoise.fasta", fasta_derep)

  # and the destination folder
  fasta_dnoise <- gsub("7_dereplication", "8_denoising", fasta_dnoise)

  # get valid FASTQ names
  to_keep <- names_checker(project_path = project_path,
                           file_folder = "7_dereplication",
                           return_fastq_names = TRUE,
                           n_reads = n_reads) %>%
    paste(., collapse = "|")

  # select cutadapt files with grepl
  fasta_derep <- fasta_derep[grepl(to_keep,  fasta_derep)] %>%
    sort()

  # select trimmed files with grepl
  fasta_dnoise <- fasta_dnoise[grepl(to_keep,  fasta_dnoise)] %>%
    sort()

  # create a read count data.frame after trimming
  read_count_dnoise <- data.frame(samples_name = gsub("_dnoise.fasta", "", basename(fasta_dnoise)),
                                  dnoise = NA)


  # create a file to store the messages from cutadapt
  file.create(file.path(project_path, "log_files", "8_dnoise.txt"))

  # for loop to store make trimming
  for(i in 1:length(fasta_dnoise)){

    if(i == 1){
      write.table(paste("###",
                        gsub("_dnoise.fasta", "", basename(fasta_dnoise[i])) ,
                        paste(rep("-", 100), collapse = "")),
                  file = file.path(project_path, "log_files", "header_primers.txt"),
                  row.names = FALSE,
                  quote = FALSE,
                  col.names = FALSE)
    } else {
      write.table(paste("\n\n\n",
                        "###",
                        gsub("_dnoise.fasta", "",basename(fasta_dnoise[i])),
                        paste(rep("-", 100), collapse = "")),
                  file = file.path(project_path, "log_files", "header_primers.txt"),
                  row.names = FALSE,
                  quote = FALSE,
                  col.names = FALSE)
    }




    file.append(file.path(project_path, "log_files", "8_dnoise.txt"),
                file.path(project_path, "log_files", "header_primers.txt"))


    cmd_den <- paste("--cluster_unoise",
                     fasta_derep[i],
                     if(isTRUE(size_in)){"--sizein"},
                     if(isTRUE(size_out)){"--sizeout"},
                     "--minsize", minsize,
                     "--unoise_alpha", unoise_alpha,
                     "--centroids",
                     fasta_dnoise[i],
                     "--fasta_width", fasta_width,
                     if(!is.null(vsearch_arguments)){vsearch_arguments},
                     sep=" ")

    # paste("-cluster_unoise ", new_names, " -centroids ", folder, "/_data/_temp_ESV.txt -unoise_alpha ", unoise_alpha, " --minsize 1 --fasta_width 0", sep="")}

    dnoise_stdout <- system2(vsearch_path,
            cmd_den,
            stdout = TRUE)



    writeLines(dnoise_stdout, file.path(project_path, "log_files", "8_dnoise_temp.txt"))

    file.append(file.path(project_path, "log_files", "8_dnoise.txt"),
                file.path(project_path, "log_files", "8_dnoise_temp.txt"))




    # populate the read counts file
    read_count_dnoise[i, 1] <- gsub("_dnoise.fasta", "", basename(fasta_dnoise[i]))
    
    bio_upload <- Biostrings::readDNAStringSet(fasta_dnoise[i])
    
    if(length(bio_upload) == 0){
      # print the progress
      message(paste(gsub("_dnoise.fasta", "", basename(fasta_dnoise[i])), ": ", i, " of ", length(fasta_dnoise), sep = ""))
      
      next()
    } else {
      read_count_dnoise[i, 2] <-  bio_upload  %>%
        names() %>%
        strsplit(., "=") %>%
        do.call(rbind, .) %>%
        as.data.frame() %>%
        dplyr::pull(2) %>%
        as.numeric() %>%
        sum()
      
      
      # print the progress
      message(paste(gsub("_dnoise.fasta", "", basename(fasta_dnoise[i])), ": ", i, " of ", length(fasta_dnoise), sep = ""))
    }
    
    
  }


  # remove temporary files to keep the log folder clean
  file.remove(file.path(project_path, "log_files", "header_primers.txt"))
  file.remove(file.path(project_path, "log_files", "8_dnoise_temp.txt"))

  # add read counts to the log file
  read_count_df <- readxl::read_excel(file.path(project_path, "log_files", "0_read_count.xlsx"))

  if("dnoise" %in% colnames(read_count_df)){
    col_to_remove <- colnames(read_count_df)[which(colnames(read_count_df) == "dnoise"):ncol(read_count_df)]
    
    # remove columns with dplyr
    read_count_df %>%
      dplyr::select(-dplyr::any_of(col_to_remove)) %>%
      dplyr::left_join(., read_count_dnoise, by = "samples_name") %>%
      writexl::write_xlsx(., file.path(project_path, "log_files", "0_read_count.xlsx"))
  } else{
    read_count_df %>%
      dplyr::left_join(., read_count_dnoise, by = "samples_name") %>%
      writexl::write_xlsx(., file.path(project_path, "log_files", "0_read_count.xlsx"))
  }
}
