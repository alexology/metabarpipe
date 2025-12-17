#' FASTQ from ENA
#'
#' @description
#' This function download
#'
#' @param x data.frame, basically the tsv downloaded from the ENA website.
#' @param project_path The path to the project folder.
#'
#' @export
#'
#' @importFrom dplyr select slice pull
#' @importFrom curl curl_fetch_disk


ena_request <- function(x = NULL,
                        project_path = NULL){

  # get the current working directory to further restore the results
  current_wd <- getwd()

  # set working directory to store the fastq
  setwd(file.path(project_path, "0_raw_files"))

  # list the files in the working directory
  fl <- list.files()


  # make a loop to download the ith fastq
  for(i in 1:nrow(x)){
    # get the accession number of the sample
    SRR <- x %>%
      dplyr::select(run_accession) %>%
      dplyr::slice(i) %>%
      dplyr::pull()

    # get ftp addresses
    fastq_ftp <- x %>%
      dplyr::select(fastq_ftp) %>%
      dplyr::slice(i) %>%
      dplyr::pull() %>%
      strsplit(., ";") %>%
      unlist()
    
    
    # request for R1
    file_r1 <- fastq_ftp[1]

    # request for R2
    file_r2 <- fastq_ftp[2]

    # check if R1 files are already present
    if(any(fl %in% basename(file_r1))){
      # print the progression
      message(paste(SRR, " already present: ", i, " of ", nrow(x),  sep = ""))
    } else {
      # make the request for R1
      curl::curl_fetch_disk(file_r1, basename(file_r1))
    }



    # check if R2 files are already present
    if(any(fl %in% basename(file_r2))){
      # print the progression
      message(paste(SRR, " already present: ", i, " of ", nrow(x),  sep = ""))
      next()
    }

    # make the request for R2
    curl::curl_fetch_disk(file_r2, basename(file_r2))

    # print the progression
    message(paste(SRR, ": ", i, " of ", nrow(x),  sep = "") )

  }

  # restore the working directory
  setwd(current_wd)
}
