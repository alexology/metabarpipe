#' OTU clustering with swarm
#'
#' @description
#' This function perfroms the OTU clustering with \code{swarm}
#'
#' @param asv_folder Path to the folder with ASV results. 
#' @param project_path Path to the project folder.
#' @param swarm_path Path to swarm.
#' @param param_d Distance parameter.
#' @param param_t Parameter t.
#' @param swarm_arguments Further arguments to be passed to \code{swarm}.
#'
#' @export
#'
#' @importFrom tibble as_tibble
#' @importFrom dplyr arrange desc mutate inner_join select
#' @importFrom Biostrings DNAStringSet width writeXStringSet
#' @importFrom readxl read_excel
#' @importFrom writexl write_xlsx
#' @importFrom utils menu read.csv
#' @importFrom rlang .data

swarm <- function(asv_folder = "13_taxonomic_assignment",
                  project_path = NULL,
                  swarm_path = NULL,
                  param_d = 13,
                  param_t = 4,
                  swarm_arguments = NULL){

  # get the asv_table from user specified folder
  file_list <- list.files(file.path(project_path, asv_folder),
                          full.names = TRUE,
                          pattern = ".xlsx")
  
  # if more then 1 file let the user choose
  if(length(file_list) > 1){
    choice_index <- menu(basename(file_list), title = "Chose the file")
    file_list <- file_list[choice_index]
  }
  
  # get the file
  asv_table <- readxl::read_excel(file_list) 
  
  # avoid R CMD checks
  asv_table$ASV_id <- paste("ASV", 1:nrow(asv_table), sep = "_")
  
  # transform the dataset to a format suitable for swarm
  otu_table_aggregated_swarm <- data.frame(ASV_id = paste("ASV", 1:nrow(asv_table), sep = "_"),
                                           sequence = asv_table$ASV,
                                           otu_abu = apply(asv_table[,sapply(asv_table, is.numeric)], 1, sum)) %>%
    tibble::as_tibble() %>%
    dplyr::arrange(dplyr::desc(.data$otu_abu))
  
  # to avoid R CMD issues
  n <- nrow(otu_table_aggregated_swarm)
  otu_table_aggregated_swarm$swarm_id <- paste("seqID", 1:n, "_", otu_table_aggregated_swarm$otu_abu, sep = "") 

  # prepare the sequences as a DNAStringSet
  otu_swarm_fasta <- Biostrings::DNAStringSet(otu_table_aggregated_swarm$sequence)
  names(otu_swarm_fasta) <- otu_table_aggregated_swarm$swarm_id

  # get the max sequence length
  max_width <- max(Biostrings::width(otu_swarm_fasta))

  # work in tempdir
  out_dir <- file.path(project_path, "swarm")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # write fasta
  Biostrings::writeXStringSet(otu_swarm_fasta,
                  filepath = file.path(out_dir,
                                       "asv_swarm.fasta"),
                  width = max_width+1)

  path_representative <- file.path(out_dir, "cluster_representatives.fasta")

  path_list <- file.path(out_dir, "cluster_list.txt")

  path_input <- file.path(out_dir, "asv_swarm.fasta")

  # run swarm
  vout <- system2(file.path(swarm_path),
                  args = c("-d", param_d,
                           "-t", param_t,
                           if(!is.null(swarm_arguments)) {swarm_arguments},
                           "-w", path_representative,
                           path_input,
                           "-o", path_list),
                  stdout = TRUE,
                  stderr = TRUE)

  if (!file.exists(path_list)) {
    stop("vsearch did not create: ", path_list,
         "\n\nvsearch output:\n", paste(vout, collapse = "\n"))
  }

  # leggo i risultati di swarm e assegno ogni ASV alla proprioa OTU
  otu_id_swarm <- read.csv(file.path(out_dir, "cluster_list.txt"), header = FALSE)
  haplo_id <- otu_id <- c()

  for(i in 1:nrow(otu_id_swarm)){
    haplo_id_temp <- strsplit(otu_id_swarm[i,], split = " ") %>%
      unlist()
    haplo_id_temp <- sub("_[^_]+$", "", haplo_id_temp)
    otu_id_temp <- rep(paste("OTU", i, sep = "_"), length(haplo_id_temp))
    haplo_id <- c(haplo_id, haplo_id_temp)
    otu_id <- c(otu_id, otu_id_temp)
  }

  # aggiorno la ASV table
  asv_table_swarm_def <-  otu_table_aggregated_swarm 
  
  # to avoid R CMD issues
  asv_table_swarm_def$swarm_id <- sub("_[^_]+$", "",  asv_table_swarm_def$swarm_id)
  
  asv_table_swarm_def <- dplyr::inner_join(data.frame(otu_id_swarm = otu_id, swarm_id = haplo_id), asv_table_swarm_def, by = "swarm_id")  %>%
    tibble::as_tibble() %>%
    dplyr::select("swarm_id", "otu_id_swarm", "ASV_id") %>%
    dplyr::inner_join(asv_table, by = "ASV_id") %>%
    dplyr::select(-c("swarm_id", "ASV_id"))

  # get the basename
  file_name <-  strsplit(basename(file_list), "\\.") %>%
    unlist()
  
  # build the file name
  file_name <- paste(file_name[1],
                     paste("otu_swarm", ".", sep = ""),
                     file_name[2],
                     sep = "")
  
  writexl::write_xlsx(x = asv_table_swarm_def,
                      path = file.path(out_dir, file_name)) 

  asv_table_swarm_def
}
