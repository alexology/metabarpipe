#' ASV table
#'
#' @description
#' This function creates an ASV table.
#'
#' @param project_path The path to the project folder.
#'
#' @export
#'
#' @importFrom dplyr pull bind_rows
#' @importFrom tibble as_tibble
#' @importFrom tidyr pivot_wider
#' @importFrom writexl write_xlsx
#' @importFrom Biostrings readDNAStringSet
#' @importFrom rlang .data

build_asv_table <- function(project_path = NULL){

  # get the file list to trim
  fasta_chimera <- list.files(file.path(project_path, "9_chimera"), full.names = TRUE) %>%
    sort()

  # site names
  site_names <- gsub("_chimera.fasta", "", basename(fasta_chimera))


  # create a list to store the results
  asv_list <- vector(mode = "list", length = length(site_names))

  # loop to retrieve the results
  for(i in 1:length(site_names)){
    # get the fasta file
    temp_fasta <- Biostrings::readDNAStringSet(filepath = fasta_chimera[i],
                                               format = "fasta")
    
    if(length(temp_fasta) == 0){
      next()
    } else {
      # get n_reads
      temp_n_reads <- strsplit(names(temp_fasta), ";") 
      
      # to avoid no visible binding for global variable '.'
      temp_n_reads <- do.call(rbind, temp_n_reads) %>%
        as.data.frame() %>%
        dplyr::pull(2) 
      
      # to avoid no visible binding for global variable '.'
      temp_n_reads <- gsub("size=", "", temp_n_reads) %>%
        as.numeric()
      
      # create the data.frame
      temp_df <- data.frame(sample = site_names[i],
                            ASV = as.character(temp_fasta),
                            n_reads = temp_n_reads) %>%
        tibble::as_tibble()
      
      
      # assign df to list
      asv_list[[i]] <- temp_df
      
    }
    
    
  }

  wide_df <- dplyr::bind_rows(asv_list) %>%
        tidyr::pivot_wider(id_cols = .data$ASV,
                names_from = .data$sample,
                values_from = .data$n_reads,
                values_fill = 0,
                values_fn = sum)

  writexl::write_xlsx(wide_df,
                      path = file.path(project_path, "10_asv_table", "asv_table.xlsx"))

  wide_df

}
