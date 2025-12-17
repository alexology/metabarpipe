#' Merge asv tables
#'
#' @description
#' This function merge multiple ASV tables
#'
#' @param ... asv_tables.
#'
#' @export
#'
#' @importFrom readxl read_excel
#' @importFrom dplyr rename filter pull
#' @importFrom tidyr pivot_longer


merge_asv_tables <-  function(...){

  # store paths as a list
  asv_table_list <- list(...)

  # get asv tables
  # asv_table_list <- lapply(tables, function(x) readxl::read_excel(file.path(x, "12_results", "asv_table_replicate_merged.xlsx")))

  # count the number of input ASV
  input_ASV <- sapply(asv_table_list, function(x) nrow(x)) %>%
    sum()

  # transform into long format
  asv_table_list <- lapply(asv_table_list,
                           function(x) tidyr::pivot_longer(data = x, cols = -ASV))


  # table to return
  to_return <- do.call(rbind, asv_table_list) %>%
    dplyr::rename(samples = name,
                  n_reads = value) %>%
    dplyr::filter(n_reads > 0)

  # count the number of ASV of the output file
  output_ASV <- to_return %>%
    dplyr::pull(ASV) %>%
    unique() %>%
    length()

  message(paste(input_ASV, "input ASV -", output_ASV, "output ASV"))

  to_return
}
