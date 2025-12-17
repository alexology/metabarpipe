#' Quantify the primer orientation
#'
#' @description
#' This function assess the primer forms. this is useful to take the reverse
#' of the complement in a later step.
#'
#'
#' @param project_path Path to the project folder.
#' @param primer_fwd Forward primer.
#' @param primer_rvr Reverse primer.
#' @param n The number of reads to subset for calculating primer occurrence.
#' Default to 1e4.
#'
#' @details
#' The function \code{primer_orientation} compare the read names of r1 and r2 on a
#' random subset of the reads. Increase \code{n} if you need more accurate results.
#'
#' @export
#'
#' @importFrom ShortRead FastqSampler yield
#' @importFrom tibble column_to_rownames rownames_to_column as_tibble
#' @importFrom writexl write_xlsx


primer_orientation <- function(project_path = NULL,
                               primer_fwd = NULL,
                               primer_rvr = NULL,
                               n = 1e4){

  # get all the primer forms that can be present in the fastq
  primer_fwd_form <- prmr_funcn(primer_fwd)
  primer_rvr_form <- prmr_funcn(primer_rvr)

  # combine the results
  primer_form <- c(primer_fwd_form, primer_rvr_form)

  # get names of r1 files
  samples_r1 <- list.files(file.path(project_path, "1_demultiplexed"),
                           pattern = "_r1.fastq",
                           full.names = TRUE) %>%
    sort()

  # get names of r2 files
  samples_r2 <- list.files(file.path(project_path, "1_demultiplexed"),
                           pattern = "_r2.fastq",
                           full.names = TRUE) %>%
    sort()

  # sample names
  sample_names <- gsub("_r1.fastq", "", basename(samples_r1))


  # create a data.frame to store the results of primer detection
  primer_count_r1 <- data.frame(primer = primer_form,
                                matrix(NA,
                                       ncol = length(samples_r1)))


  # set names for the log files
  names(primer_count_r1)[-1] <- sample_names


  # the same for r2
  primer_count_r2 <- primer_count_r1

  # count the number of times a primer appears in a subsample of the fastqs
  # this is necessary to check if there is the need to revcomp the fastq later
  for(i in 1:length(sample_names)){

    # set the subsampling of the fastq to speed up the process
    temp_r1 <- FastqSampler(samples_r1[i], n = 1e4)
    temp_r2 <- FastqSampler(samples_r2[i], n = 1e4)

    # get the sample
    temp_r1_yield <- yield(temp_r1)
    temp_r2_yield <- yield(temp_r2)

    # count r1 and r2 primers
    primer_count_r1[,i+1] <- sapply(primer_form,
                                    vcnt_funcn,
                                    filt_seq = as.character(temp_r1_yield@sread))
    primer_count_r2[,i+1] <- sapply(primer_form,
                                    vcnt_funcn,
                                    filt_seq = as.character(temp_r2_yield@sread))

    # close connections
    close(temp_r1)
    close(temp_r2)

    # remove the object from the R workspace
    rm(temp_r1)
    rm(temp_r2)

    # print the progress
    print(paste(sample_names[i], ": ", i, " of ", length(sample_names), sep = ""))
  }


  # save the results in the log file
  primer_count_r1 %>%
    tibble::column_to_rownames("primer") %>%
    t() %>%
    as.data.frame() %>%
    tibble::rownames_to_column("sample") %>%
    tibble::as_tibble()  %>%
    writexl::write_xlsx(x = .,
                        file.path(project_path,
                                  "log_files",
                                  "2_primer_count_original_fastq_R1.xlsx"))


  primer_count_r2 %>%
    tibble::column_to_rownames("primer") %>%
    t() %>%
    as.data.frame() %>%
    tibble::rownames_to_column("sample") %>%
    tibble::as_tibble()  %>%
    writexl::write_xlsx(x = .,
                        file.path(project_path,
                                  "log_files", "2_primer_count_original_fastq_R2.xlsx"))


}
