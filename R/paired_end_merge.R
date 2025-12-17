#' Paired end merge
#'
#' @description
#' This function performs paired end merging.
#'
#' @param project_path The path to the project folder.
#' @param vsearch_path The path to \code{vsearch} folder.
#' @param fastq_maxdiffs From \code{vsearch}: "... specify the maximum number of non-matching nucleotides allowed in the overlap region. The default value is 10."
#' @param fastq_maxdiffpct From \code{vsearch}: "..specify the maximum percentage of non-matching nucleotides allowed in the overlap region. The default value is 100.0 \%."
#' @param fastq_minovlen From \code{vsearch}: "...specify the minimum overlap between the merged reads. The default is 10".
#' @param mergestagger From \code{vsearch}: "...allow to merge staggered read pairs. Staggered pairs
#' are pairs where the 3' end of the reverse read has an overhang to the left of the 5' end
#'  of the forward read."
#' @param vsearch_arguments Further arguments to be passed to \code{vsearch}.
#'
#'
#' @export
#'
#' @importFrom ShortRead countFastq
#' @importFrom dplyr left_join
#' @importFrom readxl read_excel
#' @importFrom writexl write_xlsx


paired_end_merge <- function(project_path = NULL,
                             vsearch_path = NULL,
                             fastq_maxdiffs = 10,
                             fastq_maxdiffpct = 100,
                             fastq_minovlen = 10,
                             mergestagger = TRUE,
                             vsearch_arguments = NULL){

  # get the path to r1 files
  r1_path <- list.files(file.path(project_path, "2_trim_by_length"),
                        full.names = TRUE,
                        pattern = "_r1_trimmed.fastq") %>%
    sort()

  # get the path to r2 files
  r2_path <- list.files(file.path(project_path, "2_trim_by_length"),
                        full.names = TRUE,
                        pattern = "_r2_trimmed.fastq") %>%
    sort()

  # select the files based on names_checker
  to_keep <- names_checker(project_path = project_path,
                           file_folder = "2_trim_by_length",
                           r1_pattern = "_r1_trimmed.fastq",
                           r2_pattern = "_r2_trimmed.fastq",
                           return_fastq_names = TRUE) %>%
    paste(., collapse = "|")



  # select r1 files with grepl
  r1_path <- r1_path[grepl(to_keep,  r1_path)] %>%
    sort()

  # select r2 files with grepl
  r2_path <- r2_path[grepl(to_keep,  r2_path)] %>%
    sort()

  # generate new names
  merged_name <- gsub("_r1_trimmed.fastq", "_merged.fastq", r1_path)

  # adjust the directory
  merged_name <- gsub("2_trim_by_length", "3_paired_end_merge", merged_name)


  # create a table to store the results
  read_count_merged <- data.frame(samples_name = NA,
                                  merged = NA)

  # create a file to store the comments from cutadapt
  file.create(file.path(project_path, "log_files", "3_paired_end_merge.txt"))

  # apply trim by length to each file using mapply
  for(i in 1:length(merged_name)){

    if(i == 1){
      write.table(paste("###",
                        basename(merged_name)[i],
                        paste(rep("-", 100), collapse = "")),
                  file = file.path(project_path, "log_files", "header_primers.txt"),
                  row.names = FALSE,
                  quote = FALSE,
                  col.names = FALSE)
    } else {
      write.table(paste("\n\n\n",
                        "###",
                        basename(merged_name)[i],
                        paste(rep("-", 100), collapse = "")),
                  file = file.path(project_path, "log_files", "header_primers.txt"),
                  row.names = FALSE,
                  quote = FALSE,
                  col.names = FALSE)
    }

    file.append(file.path(project_path, "log_files", "3_paired_end_merge.txt"),
                file.path(project_path, "log_files", "header_primers.txt"))

    # trim fastq using usearch


    cmd <- paste(" -fastq_mergepairs ",
                 r1_path[i],
                 " -reverse ",
                 r2_path[i],
                 " -fastqout ",
                 merged_name[i],
                 " -fastq_maxdiffs ",
                 fastq_maxdiffs,
                 # usearch fastq_pctid
                 " -fastq_maxdiffpct ",
                 fastq_maxdiffpct,
                 " -fastq_minovlen ",
                 fastq_minovlen,
                 # usearch fastq_nostagger
                 if(isTRUE(mergestagger)){" --fastq_allowmergestagger "},
                 if(!is.null(vsearch_arguments)){paste(" ", vsearch_arguments, " ", sep = "")},
                 sep="")


    to_write <- system2(vsearch_path,
                        cmd,
                        stdout = TRUE)



    writeLines(to_write, file.path(project_path, "log_files", "3_paired_end_merge_temp.txt"))


    # populate the read counts file
    read_count_merged[i, 1] <- gsub("_merged.fastq", "", basename(merged_name[i]))
    read_count_merged[i, 2] <- ShortRead::countFastq(merged_name[i])[1]

    file.append(file.path(project_path, "log_files", "3_paired_end_merge.txt"),
                file.path(project_path, "log_files", "3_paired_end_merge_temp.txt"))


    message(paste(gsub("_merged.fastq", "", basename(merged_name[i])), ": ", i, " of ", length(merged_name), sep = ""))
  }

  # remove temporary files to keep the log folder clean
  file.remove(file.path(project_path, "log_files", "header_primers.txt"))
  file.remove(file.path(project_path, "log_files", "3_paired_end_merge_temp.txt"))

  # add read counts to the log file
  read_count_df <- readxl::read_excel(file.path(project_path, "log_files", "0_read_count.xlsx"))

  if("merged" %in% colnames(read_count_df)){
    col_to_remove <- colnames(read_count_df)[which(colnames(read_count_df) == "merged"):ncol(read_count_df)]
    
    # remove columns with dplyr
    read_count_df %>%
      dplyr::select(-dplyr::any_of(col_to_remove)) %>%
      dplyr::left_join(., read_count_merged, by = "samples_name") %>%
      writexl::write_xlsx(., file.path(project_path, "log_files", "0_read_count.xlsx"))
  } else{
    read_count_df %>%
      dplyr::left_join(., read_count_merged, by = "samples_name") %>%
      writexl::write_xlsx(., file.path(project_path, "log_files", "0_read_count.xlsx"))
  }
  
}
