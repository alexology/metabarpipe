#' Translation filter
#'
#' @description
#' This function select sequences that do not have a stop codon.
#'
#' @param project_path The path to the project folder. The folder must already exists.
#' @param m Minimum marker length
#' @param M Maximum marker length
#' @param genetic_code Genetic code for sequences translation.
#' @param check_reverse Also check reverse-complement reading frames? Default to \code{FALSE}
#'
#' @export
#'
#' @importFrom readxl read_excel
#' @importFrom dplyr left_join pull select rename
#' @importFrom Biostrings DNAStringSet translate subseq getGeneticCode reverseComplement vcountPattern
#' @importFrom tibble as_tibble rownames_to_column
#' @importFrom rlang .data


translate_sequences <- function(project_path = NULL,
                                m = NULL,
                                M = NULL,
                                genetic_code = NULL,
                                check_reverse = FALSE){
  
  asv_table <- readxl::read_excel(file.path(project_path, "10_asv_table", "asv_table.xlsx"))
  
  sequences <- Biostrings::DNAStringSet(asv_table %>% dplyr::pull("ASV"))
  gc <- Biostrings::getGeneticCode(as.character(genetic_code), full.search = FALSE)
  
  len_ok <- Biostrings::width(sequences) >= m & Biostrings::width(sequences) <= M
  
  f1 <- suppressWarnings(Biostrings::translate(Biostrings::subseq(sequences, 1), no.init.codon = TRUE, genetic.code = gc))
  f2 <- suppressWarnings(Biostrings::translate(Biostrings::subseq(sequences, 2), no.init.codon = TRUE, genetic.code = gc))
  f3 <- suppressWarnings(Biostrings::translate(Biostrings::subseq(sequences, 3), no.init.codon = TRUE, genetic.code = gc))
  
  ok1 <- (Biostrings::vcountPattern("*", f1) == 0) & len_ok
  ok2 <- (Biostrings::vcountPattern("*", f2) == 0) & len_ok
  ok3 <- (Biostrings::vcountPattern("*", f3) == 0) & len_ok
  
  if (isTRUE(check_reverse)) {
    rc <- Biostrings::reverseComplement(sequences)
    
    r1 <- suppressWarnings(Biostrings::translate(Biostrings::subseq(rc, 1), no.init.codon = TRUE, genetic.code = gc))
    r2 <- suppressWarnings(Biostrings::translate(Biostrings::subseq(rc, 2), no.init.codon = TRUE, genetic.code = gc))
    r3 <- suppressWarnings(Biostrings::translate(Biostrings::subseq(rc, 3), no.init.codon = TRUE, genetic.code = gc))
    
    okr1 <- (Biostrings::vcountPattern("*", r1) == 0) & len_ok
    okr2 <- (Biostrings::vcountPattern("*", r2) == 0) & len_ok
    okr3 <- (Biostrings::vcountPattern("*", r3) == 0) & len_ok
    
    keep <- rowSums(data.frame(ok1, ok2, ok3, okr1, okr2, okr3)) > 0
  } else {
    keep <- rowSums(data.frame(ok1, ok2, ok3)) > 0
  }
  
  asv_table_translated <- asv_table[keep, ]
  
  writexl::write_xlsx(
    asv_table_translated,
    path = file.path(project_path, "11_translation", "asv_table_translated.xlsx")
  )
  
  asv_table_translated_count <- readxl::read_excel(
    file.path(project_path, "11_translation", "asv_table_translated.xlsx")
  ) %>%
    dplyr::select(-.data$ASV) %>%
    as.data.frame() %>%
    (\(x) colSums(as.matrix(x)))() %>%
    as.data.frame() %>%
    tibble::rownames_to_column("samples_name") %>%
    dplyr::rename(translated_count = 2) %>%
    tibble::as_tibble()
  
  read_count_df <- readxl::read_excel(file.path(project_path, "log_files", "0_read_count.xlsx"))
  
  if ("translated_count" %in% colnames(read_count_df)) {
    read_count_df[, "translated_count"] <- asv_table_translated_count[, "translated_count"]
    writexl::write_xlsx(read_count_df, file.path(project_path, "log_files", "0_read_count.xlsx"))
  } else {
    read_count_df %>%
      dplyr::left_join(asv_table_translated_count, by = "samples_name") %>%
      writexl::write_xlsx(file.path(project_path, "log_files", "0_read_count.xlsx"))
  }
  
  asv_table_translated
}
