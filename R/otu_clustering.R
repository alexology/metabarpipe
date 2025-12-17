#' OTU clustering with vsearch
#'
#' @description
#' This function perferms the OTU clustering with \code{vsearch}.
#'
#' @param asv_table ASV table.
#' @param vsearch_path Path to vsearch.
#' @param id From \code{vsearch}: "The pairwise identity is defined as the number of (matching columns) / (alignment length - terminal gaps)".
#' @param vsearch_arguments Further arguments to be passed to \code{vsearch}.
#'
#' @export
#'
#' @importFrom tibble as_tibble
#' @importFrom dplyr arrange mutate inner_join select relocate rename
#' @importFrom Biostrings DNAStringSet width writeXStringSet readDNAStringSet
#'

otu_clustering <- function(asv_table = NULL,
                           vsearch_path = NULL,
                           id = 0.97,
                           vsearch_arguments = NULL){

  asv_table <- asv_table %>%
    dplyr::mutate(ASV_id = paste("ASV", 1:nrow(asv_table), sep = "_"))

  # transform the dataset to a format suitable for vsearch
  otu_table_aggregated_vsearch <- data.frame(ASV_id = paste("ASV", 1:nrow(asv_table), sep = "_"),
                                           sequence = asv_table$ASV,
                                           otu_abu = apply(asv_table[,sapply(asv_table, is.numeric)], 1, sum)) %>%
    tibble::as_tibble() %>%
    dplyr::arrange(desc(otu_abu)) %>%
    dplyr::mutate(vsearch_id = paste("seqID", 1:nrow(.), "_", otu_abu, sep = ""))

  # prepare the sequences as a DNAStringSet
  otu_vsearch_fasta <- Biostrings::DNAStringSet(otu_table_aggregated_vsearch$sequence)
  names(otu_vsearch_fasta) <- otu_table_aggregated_vsearch$vsearch_id

  # get the max sequence length
  max_width <- max(Biostrings::width(otu_vsearch_fasta))

  # create a directory to store vsearch results
  dir.create(file.path(getwd(), "otu_clustering"))

  # write fasta
  Biostrings::writeXStringSet(otu_vsearch_fasta,
                              filepath = file.path(getwd(),
                                                   "otu_clustering",
                                                   "asv_vsearch.fasta"),
                              width = max_width+1)

  path_representative <- file.path(getwd(), "otu_clustering", "cluster_representatives.fasta") %>%
    shQuote()

  path_list <- file.path(getwd(), "otu_clustering", "cluster_list.txt") %>%
    shQuote()

  path_input <- file.path(getwd(), "otu_clustering", "asv_vsearch.fasta") %>%
    shQuote()

  # run vsearch
  system2(vsearch_path,
          args = c("--cluster_size",
                   path_input,
                   "--centroids",
                   path_representative,
                   "--id",
                   id,
                   "--uc",
                   path_list,
                   if(!is.null(vsearch_arguments)){vsearch_arguments}))


  # leggo i risultati di vsearch e assegno ogni ASV alla propria OTU
  otu_id_vsearch <- read.table(file.path(getwd(), "otu_clustering", "cluster_list.txt"), header = FALSE) %>%
    dplyr::select(V1, V9, V10) %>%
    dplyr::filter(V1 %in% c("S", "H"))

  otu_id_vsearch_centroids <- otu_id_vsearch %>%
    dplyr::filter(V1 %in% "S") %>%
    dplyr::select(-c("V1", "V10")) %>%
    dplyr::mutate(OTU = paste("OTU", 1:nrow(.), sep = "_")) %>%
    dplyr::rename(V10 = V9)


  otu_id_vsearch_hits <- otu_id_vsearch %>%
    dplyr::filter(V1 %in% "H") %>%
    dplyr::select(-c("V1")) %>%
    dplyr::inner_join(otu_id_vsearch_centroids, ., by = "V10") %>%
    dplyr::select(-V10) %>%
    dplyr::rename(V10 = V9) %>%
    dplyr::relocate(OTU) %>%
    dplyr::bind_rows(otu_id_vsearch_centroids) %>%
    dplyr::rename("vsearch_id" = V10)

  asv_table_vsearch_def <- asv_table %>%
    dplyr::inner_join(.,   otu_table_aggregated_vsearch[, c("ASV_id", "vsearch_id")], by = "ASV_id") %>%
    dplyr::select(-ASV_id) %>%
    dplyr::inner_join(., otu_id_vsearch_hits, by = "vsearch_id") %>%
    dplyr::select(-vsearch_id) %>%
    dplyr::relocate(OTU) %>%
    dplyr::arrange(OTU)

  asv_table_vsearch_def


}
