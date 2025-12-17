#' asv to fasta
#'
#' @description
#' This function trasform an ASV table in a fasta file.
#' 
#'
#' @param asv_table An asv table with a column called ASV.
#' @param folder_path The path to the folder to store the results
#' @param file_name Name of the file, without file extension. Default to \code{asv}.
#'
#'
#' @export
#'
#'
#'
#' @importFrom Biostrings DNAStringSet writeXStringSet
#'
#'

asv_table_to_fasta <- function(asv_table,
                               folder_path = ".",
                               file_name = NULL){
  
  # the following chunk of code is to create alphanumeric strings in
  # qiime2 style
  # create a set of letters and numbers
  random_letters <- c(letters, 0:9)
  
  # count the number of sequences in the asv table
  nrow_asv <- nrow(asv_table)
  
  # set the counter for the while function
  # the number of unique alphanumeric strings need to be the same
  # of the number of sequences
  # to avoid situations in which two or more random generated strings
  # are equal a while loop is created to ensure the equality
  n_hash <- 0
  
  # while loop
  while(n_hash != nrow_asv){
    # generate alphanumeric strings with length 30
    get_hash <- lapply(rep(30, nrow_asv), function(x) sample(random_letters,
                                                        size = x,
                                                        replace = TRUE,
                                                        prob = c(rep(1/26, 26), rep(1/10, 10)))) 
    
    # collapse alpha numeric letters to string, one for each sequence
    get_hash <- sapply(get_hash, function(x) paste(x, collapse = ""))
    
    # count the number of unique ash that must be equal to the number of 
    # sequences
    n_hash <- length(unique(get_hash))
  }

  
  # create the Biostring object to save later as fasta
  fasta_asv <- Biostrings::DNAStringSet(asv_table$ASV)
  
  # name the ASVs with the alphanumeric strings created above
  names(fasta_asv) <- get_hash
  
  # set the file name to asv if none is provided by the user
  file_name <- ifelse(is.null(file_name),
                      "asv.fasta",
                      paste(file_name, ".fasta", sep = ""))
  
  # remove any white space
  file_name <- gsub(" ", "_", file_name)
  
  # set the file name for saving the fasta file
  file_path <- file.path(folder_path, file_name)
  
  # save the fasta file with writeXStringSet of the Biostrings package 
  Biostrings::writeXStringSet(fasta_asv, filepath = file_path)
  
}