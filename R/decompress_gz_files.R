#' Decompress gz file
#'
#' @description
#' Decompress gz files into fastq.
#'
#' @param project_path Path to the project folder.
#'
#' @export
#'
#' @importFrom R.utils gunzip

decompress_gz_file <- function(project_path = NULL){

  # list input files
  file_list_input <- list.files(file.path(project_path, "0_raw_files"), full.names = TRUE)

  # change the input file name to fit with the output names
  file_list_output <- gsub("0_raw_files", "1_demultiplexed", file_list_input)
  file_list_output <- gsub(".gz", "", file_list_output)

  # loop for gunzip to work
  for(i in 1:length(file_list_input)){
    R.utils::gunzip(filename = file_list_input[i],
                    destname = file_list_output[i],
                    overwrite = TRUE,
                    remove = FALSE)
    message(paste(basename(file_list_input[i]), "done:", i, "of", length(file_list_input)))
  }

}
