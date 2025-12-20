#' Trim primers
#'
#' @description
#' This function trims primers from merged reads.
#'
#' @param project_path The path to the project folder.
#' @param cutadapt_path The path to \code{cutadapt} folder.
#' @param primer_fwd Forward primer.
#' @param primer_rvr Reverse primer.
#' @param m From \code{cutadapt}: "... discard processed reads that are shorter than LENGTH".
#' @param n_cores The number of cores.
#' @param n_times From \code{cutadapt}: "... number of times an adapter is removed on the same sequence".
#' @param n_errors The number of allowed errors.
#' @param discard_untrimmed From \code{cutadapt}: "... --pair-filter=any with --discard-untrimmed, the pair is discarded if one of the reads does not contain an adapter".
#' @param cutadapt_arguments Further arguments to be passed to \code{cutadapt}.
#'
#' @details
#' See cutadapt website \url{https://cutadapt.readthedocs.io/en/stable/guide.html}.
#'
#'
#' @export
#'
#' @importFrom ShortRead countFastq
#' @importFrom dplyr filter select left_join
#' @importFrom Biostrings reverseComplement DNAString
#' @importFrom readxl read_excel
#' @importFrom writexl write_xlsx


trim_primers <- function(project_path = NULL,
                         cutadapt_path = NULL,
                         primer_fwd = NULL,
                         primer_rvr = NULL,
                         m = 0,
                         n_cores = 1,
                         n_times = 1,
                         n_errors = 0,
                         discard_untrimmed = TRUE,
                         cutadapt_arguments = NULL){

  # get merged results
  merged_name <- list.files(file.path(project_path, "3_paired_end_merge"), full.names = TRUE) %>%
    sort()

  # sample names to save cutadapt results
  fastq_cutadapt <- gsub("_merged.fastq", "_merged_cutadapt.fastq", merged_name)
  fastq_cutadapt <- gsub("3_paired_end_merge", "4_cutadapt_trimmed", fastq_cutadapt)


  # create a read count data.frame after trimming
  read_count_cutadapt <- data.frame(samples_name = gsub("_merged_cutadapt.fastq", "", basename(fastq_cutadapt)),
                                    cutadapt = NA)

  # create a file to store the comments from cutadapt
  file.create(file.path(project_path, "log_files", "4_cutadapt_output.txt"))

  # remove primers with cutadapt
  for(i in 1:length(fastq_cutadapt)){

    # get r1 primer
    temp_sample_r1 <- gsub("_merged_cutadapt.fastq", "", basename(fastq_cutadapt)[i])
    temp_primers_r1 <- readxl::read_excel(file.path(project_path, "log_files", "2_primer_count_original_fastq_R1.xlsx")) %>%
      dplyr::filter(sample == temp_sample_r1) %>%
      as.data.frame() %>%
      dplyr::select(-1) %>%
      t()

    # get the primer with more reads
    temp_primers_r1 <- temp_primers_r1[which.max(temp_primers_r1[,1]),] %>%
      names()


    if(i == 1){
      write.table(paste("###",
                        gsub("_merged_cutadapt.fastq", "", basename(fastq_cutadapt[i])) ,
                        paste(rep("-", 100), collapse = "")),
                  file = file.path(project_path, "log_files", "header_primers.txt"),
                  row.names = FALSE,
                  quote = FALSE,
                  col.names = FALSE)
    } else {
      write.table(paste("\n\n\n",
                        "###",
                        gsub("_merged_cutadapt.fastq", "",basename(fastq_cutadapt[i])),
                        paste(rep("-", 100), collapse = "")),
                  file = file.path(project_path, "log_files", "header_primers.txt"),
                  row.names = FALSE,
                  quote = FALSE,
                  col.names = FALSE)
    }

    file.append(file.path(project_path, "log_files", "4_cutadapt_output.txt"),
                file.path(project_path, "log_files", "header_primers.txt"))


    if(temp_primers_r1 == primer_fwd){
      fwd_primer <- primer_fwd
      rvr_primer <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(primer_rvr)))
    } else{
      fwd_primer <- primer_rvr
      rvr_primer <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(primer_fwd)))
    }

    # set the adapter argument for cutadapt
    primers <- paste(" -g ", fwd_primer, "...", rvr_primer, " ", sep = "")


    # set the command
    # n = number of times an adapter is removed on the same sequence
    # m = discard processed reads that are shorter than LENGTH
    # e = number of errors
    # --pair-filter=any with --discard-untrimmed, the pair is discarded if
    # one of the reads does not contain an adapter…

    trim_cmd <- paste("-n", n_times,
                      "-m", m,
                      "-e", n_errors,
                      "-j", n_cores,
                      primers,
                      if(isTRUE(discard_untrimmed)){"--discard-untrimmed"},
                      if(!is.null(cutadapt_arguments)){cutadapt_arguments},
                      "-o ", fastq_cutadapt[i],
                      merged_name[i])

    to_write <- system2(cutadapt_path,
                        args = trim_cmd,
                        stdout = TRUE)

    writeLines(to_write, file.path(project_path, "log_files", "4_cutadapt_output_temp.txt"))


    # populate the read counts file
    read_count_cutadapt[i, 1] <- gsub("_merged_cutadapt.fastq", "", basename(fastq_cutadapt[i]))
    read_count_cutadapt[i, 2] <- ShortRead::countFastq(fastq_cutadapt[i])[1]

    file.append(file.path(project_path, "log_files", "4_cutadapt_output.txt"),
                file.path(project_path, "log_files", "4_cutadapt_output_temp.txt"))

    message(paste(gsub("_merged.fastq", "", basename(fastq_cutadapt[i])), ": ", i, " of ", length(fastq_cutadapt), sep = ""))
  }

  # remove temporary files to keep the log folder clean
  file.remove(file.path(project_path, "log_files", "header_primers.txt"))
  file.remove(file.path(project_path, "log_files", "4_cutadapt_output_temp.txt"))

  # add read counts to the log file
  read_count_df <- readxl::read_excel(file.path(project_path, "log_files", "0_read_count.xlsx"))

  if("cutadapt" %in% colnames(read_count_df)){
    read_count_df[, "cutadapt"] <- read_count_cutadapt[,"cutadapt"]
    writexl::write_xlsx(read_count_df, file.path(project_path, "log_files", "0_read_count.xlsx"))
  } else{
    read_count_df  %>%
      dplyr::left_join(read_count_cutadapt, by = "samples_name") %>%
      writexl::write_xlsx(file.path(project_path, "log_files", "0_read_count.xlsx"))
  }
}
