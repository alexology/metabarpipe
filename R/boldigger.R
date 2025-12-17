#' boldigger
#' 
#' @description
#' This is a wrapper of the python program boldigger3.
#' 
#' @param asv_folder Path to the folder with ASV results. 
#' @param project_path The path to the project folder.
#' @param database_nr Is a number between 1 and 7 corresponding to the seven databases BOLD v5 currently offers.
#' @param operating_mode Is a number between 1 and the corresponding to the 3 operating modes BOLD v5 currently offers.
#' @param thresholds Threshold for taxonomic identification (defaul Species: 97, Genus: 95, Family: 90, Order: 85)
#'
#' @export
#' 
#' @importFrom Biostrings readDNAStringSet writeXStringSet
#' @importFrom readxl read_excel
#' 

boldigger <- function(asv_folder = NULL,
                      project_path = NULL,
                      database_nr = 1,
                      operating_mode = 1,
                      thresholds = NULL){
  
  
  
  if(is.null(project_path)){
    project_path <- getwd()
  }
  
  # remove existing files, otherwise unexpected beahviours
  file_to_remove <- list.files(file.path(project_path, "13_taxonomic_assignment"),
                               full.names = TRUE)
  
  # remove files if any
  if(length(file_to_remove) > 0){
    file.remove(file_to_remove)
  }
  
  # on exit set the working directory on folder path
  on.exit(setwd(project_path))
  
  # list asv file
  asv_file <- list.files(file.path(project_path, asv_folder), full.names = TRUE)
  
  # read the asv file
  asv_file <- readxl::read_excel(asv_file)
  
  # create a Biostring obiect
  asv_bio <- asv_file %>%
    dplyr::pull(ASV) %>%
    Biostrings::DNAStringSet()
  
  # give names to asv
  names(asv_bio) <- paste("asv", 1:length(asv_bio), sep = "")
  
  # set the path to save fasta for boldigger
  boldigger_fasta <- file.path(project_path, "13_taxonomic_assignment", "boldigger.fasta")
  
  # write fasta to disk
  Biostrings::writeXStringSet(x = asv_bio,
                              filepath = boldigger_fasta ,
                              width = max(Biostrings::width(asv_bio)) + 1)
  
  cmd <- paste("boldigger3 identify",
               shQuote(boldigger_fasta),
               "--db", database_nr,
               "--mode", operating_mode,
               ifelse(is.null(thresholds), "", paste("--thresholds", thresholds)),
               sep = " ")
  
  system(cmd)
}
