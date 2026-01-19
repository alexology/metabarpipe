#' OTU clustering with vsearch
#'
#' @description
#' This function perferms the OTU clustering with \code{vsearch}.
#'
#' @param asv_folder Path to the folder with ASV results. 
#' @param project_path Path to the project folder.
#' @param vsearch_path Path to vsearch.
#' @param id From \code{vsearch}: "The pairwise identity is defined as the number of (matching columns) / (alignment length - terminal gaps)".
#' @param vsearch_arguments Further arguments to be passed to \code{vsearch}.
#'
#' @export
#'
#' @importFrom tibble as_tibble
#' @importFrom dplyr arrange mutate inner_join select relocate rename
#' @importFrom Biostrings DNAStringSet width writeXStringSet readDNAStringSet
#' @importFrom rlang .data
#' @importFrom readxl read_excel
#' @importFrom writexl write_xlsx
#' @importFrom utils menu

otu_clustering <- function(asv_folder = "13_taxonomic_assignment",
                           project_path = NULL,
                           vsearch_path = NULL,
                           id = 0.97,
                           vsearch_arguments = NULL){
  
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

  # transform the dataset to a format suitable for vsearch
  otu_table_aggregated_vsearch <- data.frame(ASV_id = paste("ASV", 1:nrow(asv_table), sep = "_"),
                                             sequence = asv_table$ASV,
                                             otu_abu = apply(asv_table[,sapply(asv_table, is.numeric)], 1, sum)) %>%
    tibble::as_tibble() %>%
    dplyr::arrange(dplyr::desc(.data$otu_abu))
  
  # to avoid R CMD checks
  n <- nrow(otu_table_aggregated_vsearch)
  otu_table_aggregated_vsearch$vsearch_id <- paste("seqID", 1:n, "_", otu_table_aggregated_vsearch$otu_abu, sep = "")
  
  # prepare the sequences as a DNAStringSet
  otu_vsearch_fasta <- Biostrings::DNAStringSet(otu_table_aggregated_vsearch$sequence)
  names(otu_vsearch_fasta) <- otu_table_aggregated_vsearch$vsearch_id
  
  # get the max sequence length
  max_width <- max(Biostrings::width(otu_vsearch_fasta))
  
  # work in tempdir
  out_dir <- file.path(project_path, "otu_clustering")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
 
  # write fasta
  Biostrings::writeXStringSet(x = otu_vsearch_fasta,
                              filepath = file.path(out_dir, "asv_vsearch.fasta"),
                              width = max_width+1)
  
  path_representative <- file.path(out_dir, "cluster_representatives.fasta")
  
  path_list <- file.path(out_dir, "cluster_list.txt") 
  
  path_input <- file.path(out_dir, "asv_vsearch.fasta")

  cmd_oc <- c("--cluster_size",
              path_input,
              "--centroids",
              path_representative,
              "--id",
              as.character(id),
              "--uc",
              path_list,
              if(!is.null(vsearch_arguments)){vsearch_arguments},
              sep = " ")
  
  # run vsearch
  vout <- system2(vsearch_path,
          args = cmd_oc,
          stdout = TRUE,
          stderr = TRUE)
  
  if (!file.exists(path_list)) {
    stop("vsearch did not create: ", path_list,
         "\n\nvsearch output:\n", paste(vout, collapse = "\n"))
  }
  

  # leggo i risultati di vsearch e assegno ogni ASV alla propria OTU
  otu_id_vsearch <- read.table(file.path(out_dir, "cluster_list.txt"), header = FALSE) %>%
    dplyr::select("V1", "V9", "V10") %>%
    dplyr::filter(.data$V1 %in% c("S", "H"))
  
  otu_id_vsearch_centroids <- otu_id_vsearch %>%
    dplyr::filter(.data$V1 %in% "S") %>%
    dplyr::select(-c("V1", "V10")) 
  
  n <- nrow(otu_id_vsearch_centroids)
  
  otu_id_vsearch_centroids <- otu_id_vsearch_centroids %>%
    dplyr::mutate(OTU = paste("OTU", 1:n, sep = "_")) 
  
  names(otu_id_vsearch_centroids)[names(otu_id_vsearch_centroids) == "V9"] <- "V10"

  otu_id_vsearch_hits <- otu_id_vsearch %>%
    dplyr::filter(.data$V1 %in% "H") %>%
    dplyr::select(-c("V1"))
  
  otu_id_vsearch_hits <- dplyr::inner_join(otu_id_vsearch_centroids, otu_id_vsearch_hits, by = "V10") %>%
    dplyr::select(-c("V10")) 
  
  names(otu_id_vsearch_hits)[names(otu_id_vsearch_hits) == "V9"] <- "V10"
  
  otu_id_vsearch_hits <- otu_id_vsearch_hits %>%
    dplyr::relocate("OTU") %>%
    dplyr::bind_rows(otu_id_vsearch_centroids) %>%
    dplyr::rename("vsearch_id" = "V10")
  
  asv_table_vsearch_def <- asv_table %>%
    dplyr::inner_join(otu_table_aggregated_vsearch[, c("ASV_id", "vsearch_id")], by = "ASV_id") %>%
    dplyr::select(-c("ASV_id")) %>%
    dplyr::inner_join(otu_id_vsearch_hits, by = "vsearch_id") %>%
    dplyr::select(-c("vsearch_id")) %>%
    dplyr::relocate("OTU") 
  
  # sort by OTU, not sure arrange is working
  asv_table_vsearch_def <- asv_table_vsearch_def[order(asv_table_vsearch_def$OTU, decreasing = FALSE), ]
  
  # get the basename
  file_name <-  strsplit(basename(file_list), "\\.") %>%
    unlist()
  
  # build the file name
  file_name <- paste(file_name[1],
                     paste("otu_vsearch_id_", id, ".", sep = ""),
                     file_name[2],
                     sep = "")
  
  # remove the point on id threshold
  file_name <- gsub("0.|.00", "", file_name)
  
  writexl::write_xlsx(x = asv_table_vsearch_def,
                      path = file.path(out_dir, file_name)) 
    
  asv_table_vsearch_def
}
