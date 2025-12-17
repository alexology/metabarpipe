#' OTU clustering with swarm
#'
#' @description
#' This function perfroms the OTU clustering with \code{swarm}
#'
#' @param asv_table ASV table.
#' @param swarm_path Path to swarm.
#' @param param_d Distance parameter.
#' @param param_t Parameter t.
#' @param swarm_arguments Further arguments to be passed to \code{swarm}.
#'
#' @export
#'
#' @importFrom tibble as_tibble
#' @importFrom dplyr arrange mutate inner_join select
#' @importFrom Biostrings DNAStringSet width writeXStringSet
#'

swarm <- function(asv_table = NULL,
                  swarm_path = NULL,
                  param_d = 13,
                  param_t = 4,
                  swarm_arguments = NULL){

  asv_table <- asv_table %>%
    dplyr::mutate(ASV_id = paste("ASV", 1:nrow(asv_table), sep = "_"))

  # transform the dataset to a format suitable for swarm
  otu_table_aggregated_swarm <- data.frame(ASV_id = paste("ASV", 1:nrow(asv_table), sep = "_"),
                                           sequence = asv_table$ASV,
                                           otu_abu = apply(asv_table[,sapply(asv_table, is.numeric)], 1, sum)) %>%
    tibble::as_tibble() %>%
    dplyr::arrange(desc(otu_abu)) %>%
    dplyr::mutate(swarm_id = paste("seqID", 1:nrow(.), "_", otu_abu, sep = ""))

  # prepare the sequences as a DNAStringSet
  otu_swarm_fasta <- Biostrings::DNAStringSet(otu_table_aggregated_swarm$sequence)
  names(otu_swarm_fasta) <- otu_table_aggregated_swarm$swarm_id

  # get the max sequence length
  max_width <- max(Biostrings::width(otu_swarm_fasta))

  # create a directory to store swarm results
  dir.create(file.path(getwd(), "swarm"))

  # write fasta
  Biostrings::writeXStringSet(otu_swarm_fasta,
                  filepath = file.path(getwd(),
                                       "swarm",
                                       "asv_swarm.fasta"),
                  width = max_width+1)

  path_representative <- file.path(getwd(), "swarm", "cluster_representatives.fasta") %>%
    shQuote()



  path_list <- file.path(getwd(), "swarm", "cluster_list.txt") %>%
    shQuote()


  path_input <- file.path(getwd(), "swarm", "asv_swarm.fasta") %>%
    shQuote()

  # run swarm
  system2(file.path(swarm_path, "bin/swarm.exe"),
          args = c("-d", param_d,
                   "-t", param_t,
                   if(!is.null(swarm_arguments)) {swarm_arguments},
                   "-w", path_representative,
                   path_input,
                   "-o", path_list))


  # leggo i risultati di swarm e assegno ogni ASV alla proprioa OTU
  otu_id_swarm <- read.csv(file.path(getwd(), "swarm", "cluster_list.txt"), header = FALSE)
  haplo_id <- otu_id <- c()

  for(i in 1:nrow(otu_id_swarm)){
    haplo_id_temp <- strsplit(otu_id_swarm[i,], split = " ") %>%
      unlist() %>%
      sub("_[^_]+$", "", .)
    otu_id_temp <- rep(paste("OTU", i, sep = "_"), length(haplo_id_temp))
    haplo_id <- c(haplo_id, haplo_id_temp)
    otu_id <- c(otu_id, otu_id_temp)
  }

  # aggiorno la ASV table
  asv_table_swarm_def <-  otu_table_aggregated_swarm %>%
    dplyr::mutate(swarm_id = sub("_[^_]+$", "", swarm_id)) %>%
    dplyr::inner_join(data.frame(otu_id_swarm = otu_id, swarm_id = haplo_id), ., by = "swarm_id")  %>%
    tibble::as_tibble() %>%
    dplyr::select(swarm_id, otu_id_swarm, ASV_id) %>%
    dplyr::inner_join(., asv_table, by = "ASV_id") %>%
    dplyr::select(-c("swarm_id", "ASV_id"))


  asv_table_swarm_def


}
