#' Fastq reverse complement
#'
#' @description
#' This function compute the reverse complement of a fastq.
#'
#' @param project_path The path to the project folder.
#' @param primer_fwd Forward primer.
#' @param vsearch_path The path to \code{vsearch} folder.
#'
#' @details
#' The function \code{reverse_complement} replace the original file with the reverse
#' complement. The result is irreversible.
#'
#' @export
#'
#' @importFrom dplyr select
#' @importFrom utils write.table
#' @importFrom readxl read_excel


reverse_complement <- function(project_path = NULL,
                               primer_fwd = NULL,
                               vsearch_path = NULL){
  # get primer counts
  primer_count_r1 <- readxl::read_excel(file.path(project_path, "log_files", "2_primer_count_original_fastq_R1.xlsx"))
  primer_count_r2 <- readxl::read_excel(file.path(project_path, "log_files", "2_primer_count_original_fastq_R2.xlsx"))

  # get the list of cutadapt samples
  fastq_cutadapt <- list.files(file.path(project_path, "4_cutadapt_trimmed"),
                               full.names = TRUE)

  # create a file to store the comments from cutadapt
  file.create(file.path(project_path, "log_files", "reverse_complement_report.txt"))

  # a for loop that detect the orientation and change automatically the orientation
  # original files are saved in the folder original_orientation_files
  for(i in 1:length(fastq_cutadapt)){
    # get r1 primer
    # first, get the sample name
    temp_sample_r1 <- gsub("_merged_cutadapt.fastq", "", basename(fastq_cutadapt)[i])

    # second get the list of primers for the target sample
    temp_primers_r1 <- readxl::read_excel(file.path(project_path, "log_files", "2_primer_count_original_fastq_R1.xlsx")) %>%
      dplyr::filter(sample == temp_sample_r1) %>%
      as.data.frame() %>%
      dplyr::select(-1) %>%
      t()

    # third, get the primer with more reads
    temp_primers_r1 <- temp_primers_r1[which.max(temp_primers_r1[,1]),] %>%
      names()


    # set the log file to keep track of the changes


    if(temp_primers_r1 == primer_fwd){
      # print the progress

      to_print <- paste(temp_sample_r1, ": ", i, " of ", length(fastq_cutadapt), sep = "")

      utils::write.table(to_print,
                         file = file.path(project_path, "log_files", "reverse_complement_report.txt"),
                         row.names = FALSE,
                         quote = FALSE,
                         col.names = FALSE,
                         append = TRUE)

      message(to_print)

    } else {
      original_file_r1 <- fastq_cutadapt[i]
      changed_name_r1 <- gsub("_merged_cutadapt.fastq", "_original.fastq", original_file_r1)

      file.rename(from = original_file_r1,
                  to = changed_name_r1)

      r1_cmd <- paste("-fastx_revcomp \"", changed_name_r1, "\" -fastqout \"", original_file_r1,  sep="")

      system2(vsearch_path,
              r1_cmd,
              stdout = FALSE)

      file.remove(changed_name_r1)

      to_print <- paste(temp_sample_r1, ": ", i, " of ", length(fastq_cutadapt), " - reversed", sep = "")

      utils::write.table(to_print,
                         file = file.path(project_path, "log_files", "reverse_complement_report.txt"),
                         row.names = FALSE,
                         quote = FALSE,
                         col.names = FALSE,
                         append = TRUE)


      # print the progress
      message(to_print)

    }
  }
}

