#' Build BOLD taxonomy
#'
#' @description
#' This function build the taxonomy of species in an ASV table using BOLD.
#'
#' @param folder_path Path to the folder where results are stored.
#' @param asv_table An ASV table containing the taxonomic information
#' @param query_size Maximum query size to BOLD.
#'
#' @export
#'
#' @importFrom dplyr filter pull select everything bind_rows distinct_all
#' @importFrom tidyr pivot_longer
#' @importFrom bold bold_tax_name bold_tax_id2
#' @importFrom tibble column_to_rownames as_tibble
#' @importFrom writexl write_xlsx


build_taxonomy <-  function(folder_path = NULL,
                            asv_table = NULL,
                            query_size = 500){

  # get unique taxa
  asv_unique <- asv_table %>%
    dplyr::filter(! is.na(Taxa)) %>%
    dplyr::pull(Taxa) %>%
    unique()

  # make chunks to avoid length request issues
  chunks <- split(1:length(asv_unique), ceiling(seq_along(1:length(asv_unique))/query_size))

  # get bold_taxonomy
  bold_taxonomy <- lapply(chunks, function(x)bold::bold_tax_name(asv_unique[x]))

  # transform the call to a data.frame
  bold_taxonomy <- do.call(rbind, bold_taxonomy)

  # get the tree
  bold_tree <- bold::bold_tax_id2(bold_taxonomy$taxid, includeTree = TRUE)

  # put the tree in a usable format
  # get unique input id
  taxid_unique <- bold_tree %>%
    dplyr::pull(input) %>%
    unique()

  # list to store the results
  store_list <- vector(mode = "list", length = length(taxid_unique))

  # iterate over unique ids
  for(i in 1:length(store_list)){
    store_list[[i]] <- bold_tree %>%
      dplyr::filter(input == taxid_unique[i]) %>%
      dplyr::select(tax_rank, taxon) %>%
      tibble::column_to_rownames("tax_rank") %>%
      t() %>%
      as.data.frame() %>%
      tibble::as_tibble()

  }

  # transform the list into a data.frame
  store_list <- dplyr::bind_rows(store_list) %>%
    dplyr::distinct_all() %>%
    dplyr::select(any_of(c("phylum", "class", "order", "family", "subfamily", "genus", "species")))

  # tranform the column names in a format suitable for biomonitoR
  colnames_store_list <- colnames(store_list)
  colnames(store_list) <- paste(toupper(substr(colnames_store_list, 1, 1)), substr(colnames_store_list, 2, nchar(colnames_store_list)), sep="")

  # taxonomic tree
  taxonomic_tree <-  tree_from_bold(as.data.frame(store_list))

  # taxonomic list
  taxonomic_list <- taxonomic_tree %>%
    dplyr::select(-Taxa) %>%
    tidyr::pivot_longer(dplyr::everything()) %>%
    dplyr::filter(value != "") %>%
    dplyr::distinct_all()

  writexl::write_xlsx(taxonomic_tree, file.path(folder_path, "taxonomic_tree.xlsx"))
  writexl::write_xlsx(taxonomic_list, file.path(folder_path, "taxonomic_list.xlsx"))

}
