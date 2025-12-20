#' Chimera removal
#'
#' @description
#' This function perfroms chimera removal.
#'
#' @param project_path The path to the project folder.
#' @param vsearch_path The path to \code{vsearch} folder.
#' @param size_in From \code{vsearch}: "In de novo mode, abundance annotations present in sequence headers are taken into account by default".
#' @param size_out From \code{vsearch}: "...add abundance annotations to fasta headers".
#' @param fasta_width From \code{vsearch}: "Fasta files produced by vsearch are wrapped (sequences are written on lines of integer nucleotides, 80 by default). Set that value to zero to eliminate the wrapping".
#' @param relabel From \code{vsearch}: "Relabel sequences using the prefix string and a ticker (1, 2, 3, etc.) to construct the new headers. Use --sizeout to conserve the abundance annotations".
#' @param vsearch_arguments Further arguments to be passed to \code{vsearch}.
#' @param n_reads The number of reads above which a fasta is kept for further processing.
#'
#' @export
#'
#' @importFrom dplyr pull left_join
#' @importFrom writexl write_xlsx
#' @importFrom Biostrings readDNAStringSet




remove_chimera <-  function(project_path,
                            vsearch_path,
                            size_in = TRUE,
                            size_out = TRUE,
                            fasta_width = 0,
                            relabel = NULL,
                            vsearch_arguments = NULL,
                            n_reads = 0){
  # get the file list to trim
  fasta_dnoise <- list.files(file.path(project_path, "8_denoising"), full.names = TRUE) %>%
    sort()

  # change the names
  fasta_chimera <- gsub("_dnoise.fasta", "_chimera.fasta", fasta_dnoise)

  # and the destination folder
  fasta_chimera <- gsub("8_denoising", "9_chimera", fasta_chimera)


  # get valid FASTQ names
  to_keep <- names_checker(project_path = project_path,
                           file_folder = "8_denoising",
                           return_fastq_names = TRUE,
                           n_reads = n_reads) %>%
    paste(collapse = "|")

  # select cutadapt files with grepl
  fasta_dnoise <- fasta_dnoise[grepl(to_keep,  fasta_dnoise)] %>%
    sort()

  # select trimmed files with grepl
  fasta_chimera <- fasta_chimera[grepl(to_keep,  fasta_chimera)] %>%
    sort()


  # create a read count data.frame after trimming
  read_count_chimera <- data.frame(samples_name = gsub("_dnoise.fasta", "", basename(fasta_chimera)),
                                   chimera = NA)

  # create a file to store the messages from cutadapt
  file.create(file.path(project_path, "log_files", "9_chimera.txt"))

  # for loop to store make trimming
  for(i in 1:length(fasta_chimera)){

    if(i == 1){
      write.table(paste("###",
                        gsub("_chimera.fasta", "", basename(fasta_chimera[i])) ,
                        paste(rep("-", 100), collapse = "")),
                  file = file.path(project_path, "log_files", "header_primers.txt"),
                  row.names = FALSE,
                  quote = FALSE,
                  col.names = FALSE)
    } else {
      write.table(paste("\n\n\n",
                        "###",
                        gsub("_chimera.fasta", "",basename(fasta_chimera[i])),
                        paste(rep("-", 100), collapse = "")),
                  file = file.path(project_path, "log_files", "header_primers.txt"),
                  row.names = FALSE,
                  quote = FALSE,
                  col.names = FALSE)
    }




    file.append(file.path(project_path, "log_files", "9_chimera.txt"),
                file.path(project_path, "log_files", "header_primers.txt"))



    cmd_chi <- paste(" --uchime3_denovo ",
                     fasta_dnoise[i],
                     if(isTRUE(size_in)){"--sizein"},
                     if(isTRUE(size_out)){"--sizeout"},
                     "--fasta_width", fasta_width,
                     "--nonchimeras ",
                     fasta_chimera[i],
                     if(!is.null(relabel)){paste("--relabel", relabel)},
                     if(!is.null(vsearch_arguments)){vsearch_arguments},
                     sep = " ")



    chimera_stdout <- system2(vsearch_path,
                              cmd_chi,
                              stdout = TRUE)



    writeLines(chimera_stdout, file.path(project_path, "log_files", "9_chimera_temp.txt"))

    file.append(file.path(project_path, "log_files", "9_chimera.txt"),
                file.path(project_path, "log_files", "9_chimera_temp.txt"))




    # populate the read counts file
    read_count_chimera[i, 1] <- gsub("_chimera.fasta", "", basename(fasta_chimera[i]))
    
    bio_upload <- Biostrings::readDNAStringSet(fasta_chimera[i])
    
    if(length(bio_upload) == 0){
      # print the progress
      message(paste(gsub("_chimera.fasta", "", basename(fasta_chimera[i])), ": ", i, " of ", length(fasta_chimera), sep = ""))
      next()
    } else {
      
      read_count_chimera_temp <- strsplit(names(bio_upload), "=")
      read_count_chimera[i, 2] <- do.call(rbind, read_count_chimera_temp) %>%
        as.data.frame() %>%
        dplyr::pull(2) %>%
        as.numeric() %>%
        sum()
      
      
      # print the progress
      message(paste(gsub("_chimera.fasta", "", basename(fasta_chimera[i])), ": ", i, " of ", length(fasta_chimera), sep = ""))
      
    }
    


  }


  # remove temporary files to keep the log folder clean
  file.remove(file.path(project_path, "log_files", "header_primers.txt"))
  file.remove(file.path(project_path, "log_files", "9_chimera_temp.txt"))

  # add read counts to the log file
  read_count_df <- readxl::read_excel(file.path(project_path, "log_files", "0_read_count.xlsx"))

  if("chimera" %in% colnames(read_count_df)){
    col_to_remove <- colnames(read_count_df)[which(colnames(read_count_df) == "chimera"):ncol(read_count_df)]
    
    # remove columns with dplyr
    read_count_df %>%
      dplyr::select(-dplyr::any_of(col_to_remove)) %>%
      dplyr::left_join(read_count_chimera, by = "samples_name") %>%
      writexl::write_xlsx(file.path(project_path, "log_files", "0_read_count.xlsx"))
  } else{
    read_count_df %>%
      dplyr::left_join(read_count_chimera, by = "samples_name") %>%
      writexl::write_xlsx(file.path(project_path, "log_files", "0_read_count.xlsx"))
  }
}
