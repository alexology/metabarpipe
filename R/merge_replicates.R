#' Merge sample replicates
#'
#' @description
#' This function merge multiple replicates of the same sample.
#'
#' @param project_path Path to the project folder.
#' @param asv_table Can be \code{chimera_removed} for the asv_table generated after
#' chimera removal or \code{chimera} for the asv_table generated after protein translation.
#' @param pattern The pattern to remove. See details.
#' @param percentage_to_remove Remove ASV with a relative abundance lower than this threshold.
#'
#' @details
#' Replicates of the same sample are usually stored with the name of the sample and
#' an indication of the replicate (e.g. replicate1, replicate2). Please provide a pattern
#' without the replicate number (e.g. \code{pattern = "replicate"}).
#'
#' @export
#'
#' @importFrom readxl read_excel
#' @importFrom tibble column_to_rownames rownames_to_column as_tibble
#' @importFrom dplyr mutate group_by
#' @importFrom writexl write_xlsx




merge_replicates <- function(project_path = NULL,
                             asv_table = "translated",
                             pattern = "_replicate",
                             percentage_to_remove = 0){

  if(identical(asv_table, "translated")){
    asv_table <- file.path(project_path, "11_translation", "asv_table_translated.xlsx") %>%
      readxl::read_excel()
  }

  if(identical(asv_table, "chimera_removed")){
    asv_table <- file.path(project_path, "10_asv_table", "asv_table.xlsx") %>%
      readxl::read_excel()
  }


  # transpose the asv_table for further calculations
  merged_table <- asv_table %>%
    as.data.frame() %>%
    tibble::column_to_rownames("ASV") %>%
    t() %>%
    as.data.frame()

  # calculate the relative abundance of each ASV
  merged_table_percentage <- merged_table %>%
    sweep(., 1, apply(., 1, sum), "/") * 100

  # set ASV with abundance lower than a threshold to 0
  merged_table_percentage[merged_table_percentage < percentage_to_remove] <- 0
  merged_table_percentage[merged_table_percentage >= percentage_to_remove] <- 1


  # remove ASV below the threshold
  merged_table <- merged_table * merged_table_percentage

  # prepare merged table for saving
  merged_table <- merged_table %>%
    tibble::rownames_to_column("samples") %>%
    tibble::as_tibble() %>%
    dplyr::mutate(samples = gsub(paste(pattern, ".", sep = ""), "", samples)) %>%
    dplyr::group_by(samples) %>%
    dplyr::summarise_all(function(...)if(all(...> 0)){sum(...)} else {0}) %>%
    as.data.frame() %>%
    tibble::column_to_rownames("samples") %>%
    t() %>%
    as.data.frame()



  # remove empty rows
  merged_table <- merged_table[rowSums(merged_table) > 0,] %>%
    tibble::rownames_to_column("ASV") %>%
    tibble::as_tibble()

  # save to disk
  writexl::write_xlsx(merged_table,
                      file.path(project_path, "12_results", "asv_table_replicate_merged.xlsx"))


}
