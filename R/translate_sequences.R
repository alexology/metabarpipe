#' Translation filter
#'
#' @description
#' This function select sequences that do not have a stop codon.
#'
#' @param project_path The path to the project folder. The folder must already exists.
#' @param m Minimum marker length
#' @param M Maximum marker length
#' @param genetic_code Gentic code for sequences translation. See
#'
#' @export
#'
#' @importFrom readxl read_excel
#' @importFrom dplyr pull select rename
#' @importFrom Biostrings DNAStringSet translate subseq getGeneticCode reverse vcountPattern
#' @importFrom tibble as_tibble rownames_to_column


translate_sequences <- function(project_path = NULL,
                                m = NULL,
                                M = NULL,
                                genetic_code = NULL){

  # read chimera removal results
  asv_table <- readxl::read_excel(file.path(project_path, "10_asv_table", "asv_table.xlsx"))

  # genetic code from number to character. This is a Biostrings requirement
  genetic_code <- as.character(genetic_code)


  # remove sequences with stop codon
  # get the sequences
  sequences <- asv_table %>%
    dplyr::pull("ASV")

  # transform sequences to DNAStringSet
  sequences_biostrings <- Biostrings::DNAStringSet(sequences)

  # select sequence length(s). since we want to work on haplotypes, we keep 421
  selected_widths <- m:M

  # set a vector to recover coding sequences
  select_false <- rep(FALSE, length(sequences))

  # for loop to get the results
  for(i in 1:length(selected_widths)){

    # translate and filter
    select_temp_1 <- Biostrings::DNAStringSet(Biostrings::subseq(sequences, 1)) %>%
      Biostrings::translate(., no.init.codon = TRUE, genetic.code = Biostrings::getGeneticCode(genetic_code, full.search = FALSE)) %>%
      suppressWarnings()
    select_temp_2 <- Biostrings::DNAStringSet(Biostrings::subseq(sequences, 2)) %>%
      Biostrings::translate(., no.init.codon = TRUE, genetic.code = Biostrings::getGeneticCode(genetic_code, full.search = FALSE)) %>%
      suppressWarnings()
    select_temp_3 <- Biostrings::DNAStringSet(Biostrings::subseq(sequences, 3)) %>%
      Biostrings::translate(., no.init.codon = TRUE, genetic.code = Biostrings::getGeneticCode(genetic_code, full.search = FALSE)) %>%
      suppressWarnings()

    select_temp_r1 <- Biostrings::DNAStringSet(Biostrings::subseq(sequences, 1)) %>%
      Biostrings::reverse() %>%
      Biostrings::translate(., no.init.codon = TRUE, genetic.code = Biostrings::getGeneticCode(genetic_code, full.search = FALSE)) %>%
      suppressWarnings()

    select_temp_r2 <- Biostrings::DNAStringSet(Biostrings::subseq(sequences, 2)) %>%
      Biostrings::reverse() %>%
      Biostrings::translate(., no.init.codon = TRUE, genetic.code = Biostrings::getGeneticCode(genetic_code, full.search = FALSE)) %>%
      suppressWarnings()

    select_temp_r3 <- Biostrings::DNAStringSet(Biostrings::subseq(sequences, 3)) %>%
      Biostrings::reverse() %>%
      Biostrings::translate(., no.init.codon = TRUE, genetic.code = Biostrings::getGeneticCode(genetic_code, full.search = FALSE)) %>%
      suppressWarnings()

    select_temp_1 <- (Biostrings::vcountPattern("*", select_temp_1)==0) & (nchar(sequences) == selected_widths[i])
    select_temp_2 <- (Biostrings::vcountPattern("*", select_temp_2)==0) & (nchar(sequences) == selected_widths[i])
    select_temp_3 <- (Biostrings::vcountPattern("*", select_temp_3)==0) & (nchar(sequences) == selected_widths[i])
    select_temp_r1 <- (Biostrings::vcountPattern("*", select_temp_r1)==0) & (nchar(sequences) == selected_widths[i])
    select_temp_r2 <- (Biostrings::vcountPattern("*", select_temp_r2)==0) & (nchar(sequences) == selected_widths[i])
    select_temp_r3 <- (Biostrings::vcountPattern("*", select_temp_r3)==0) & (nchar(sequences) == selected_widths[i])

    select_temp <- data.frame(select_temp_1,
                              select_temp_2,
                              select_temp_3,
                              select_temp_r1,
                              select_temp_r2,
                              select_temp_r3) %>%
      apply(., 1, sum)

    # remove unwanted ASV
    select_temp <- select_temp > 0

    # update select false
    select_false <- select_temp | select_false

  }

  # remove non coding sequences
  asv_table_translated <- asv_table[select_false, ]


  # save to disk as an excel
  writexl::write_xlsx(asv_table_translated,
                      path = file.path(project_path, "11_translation", "asv_table_translated.xlsx"))

  # count the number of reads
  asv_table_translated_count <- file.path(project_path, "11_translation", "asv_table_translated.xlsx") %>%
    readxl::read_excel() %>%
    dplyr::select(-ASV) %>%
    apply(., 2, sum) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("samples_name") %>%
    dplyr::rename(translated_count = 2) %>%
    tibble::as_tibble()

  # add read counts to the log file
  read_count_df <- readxl::read_excel(file.path(project_path, "log_files", "0_read_count.xlsx"))

  if("translated_count" %in% colnames(read_count_df)){
    read_count_df[, "translated_count"] <- asv_table_translated_count[,"translated_count"]
    writexl::write_xlsx(read_count_df, file.path(project_path, "log_files", "0_read_count.xlsx"))
  } else{
    read_count_df  %>%
      dplyr::left_join(., asv_table_translated_count, by = "samples_name") %>%
      writexl::write_xlsx(., file.path(project_path, "log_files", "0_read_count.xlsx"))
  }

  asv_table_translated

}
