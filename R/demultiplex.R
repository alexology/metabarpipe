#' Demultiplex multiplexed fastq
#'
#' @description
#' This function demultiplex multiplexed data-
#'
#' @param project_path The path to the project folder.
#' @param r1 Path to the multiplexed r1 file.
#' @param r2 Path to the multiplexed r2 file.
#' @param Indices combinations.
#' @param Indices used.
#' @param Output folder.
#'
#' @export
#'
#'
#'
#' @importFrom ShortRead FastqStreamer yield
#' @importFrom Biostrings width writeXStringSet DNAString


demultiplex <- function(project_path = NULL,
                        r1 = NULL,
                        r2 = NULL,
                        combinations = NULL,
                        indices = NULL,
                        folder = NULL){

  # set the function to import chunks of fastq in R
  strm_r1 <- ShortRead::FastqStreamer(r1)
  strm_r2 <- ShortRead::FastqStreamer(r2)

  # set the counter to get information about the chunk
  counter <- 0

  # create names for storing the results
  r1_names <- apply(combinations, 1, function(x) file.path(project_path, folder, paste("s", x[1], x[2], "r1.fastq", sep = "_")) )
  r2_names <- apply(combinations, 1, function(x) file.path(project_path, folder, paste("s", x[1], x[2], "r2.fastq", sep = "_")) )

  # create emtpy files where store the results
  cat("creating empty files...\n")
  invisible(sapply(r1_names, file.create))
  invisible(sapply(r2_names, file.create))

  repeat {

    # load r1 and r2
    fq_1 <- ShortRead::yield(strm_r1)

    if(length(fq_1) == 0) {break()}

    ## process chunk
    fq_2 <- ShortRead::yield(strm_r2)

    # check that sequences are longer than 6,  just to avoid issue with
    # VcountPattern
    r1_int <- Biostrings::width(fq_1@sread) > 6
    r2_int <- Biostrings::width(fq_2@sread) > 6

    # get the subsets with sequences greater than 6
    r1_sub <- Biostrings::subseq(fq_1@sread[which((r1_int*r2_int) > 0)], 1, 7)
    r2_sub <- Biostrings::subseq(fq_2@sread[which((r1_int*r2_int) > 0)], 1, 7)

    counter <- counter + 1
    message(paste("chunk", counter, "started..."))

    for(i in 1:nrow(combinations)){
      fwd <- combinations[i, 1]
      rvr <- combinations[i, 2]

      delete_1 <- indices[indices$ID %in% fwd, 2]
      delete_2 <- indices[indices$ID %in% rvr, 2]

      # index_1
      temp_i1 <- Biostrings::DNAString(indices[indices$ID %in% fwd, 1])

      # use DECIPHER to disambiguate
      temp_i1 <- DECIPHER::Disambiguate(Biostrings::DNAStringSet(temp_i1)) %>%
        unlist()
      
      # index_2
      temp_i2 <- Biostrings::DNAString(indices[indices$ID %in% rvr, 1])

      # use DECIPHER to disambiguate
      temp_i2 <- DECIPHER::Disambiguate(Biostrings::DNAStringSet(temp_i2)) %>%
        unlist()
      
      # get the reads with the specified tag
      temp_count_1 <- lapply(1:length(temp_i1), function(x) Biostrings::vcountPattern(temp_i1[[x]],
                                                                                      Biostrings::subseq(r1_sub, 1, length(temp_i1[[x]])),
                                                                                      fixed = TRUE))
      temp_count_1 <- Reduce("+", temp_count_1)
      

      temp_count_2 <- lapply(1:length(temp_i2), function(x) Biostrings::vcountPattern(temp_i2[[x]],
                                                                                      Biostrings::subseq(r2_sub, 1, length(temp_i2[[x]])),
                                                                                      fixed = TRUE))
      temp_count_2 <- Reduce("+", temp_count_2)
      

      temp_r1 <- fq_1@sread[which((r1_int*r2_int) > 0)][which(temp_count_1*temp_count_2 > 0)]
      
      if(length(temp_r1) < 1){
        next()
      }
      
      temp_r2 <- fq_2@sread[which((r1_int*r2_int) > 0)][which(temp_count_1*temp_count_2 > 0)]

      temp_r1 <- Biostrings::subseq(temp_r1, delete_1 + 1)
      temp_r2 <- Biostrings::subseq(temp_r2, delete_2 + 1)

      # r1_fq <- c(r1_fq, temp_r1)
      # r2_fq <- c(r2_fq, temp_r2)

      temp_id_r1 <- fq_1@id[which((r1_int*r2_int) > 0)][which(temp_count_1*temp_count_2 > 0)]
      temp_id_r2 <- fq_2@id[which((r1_int*r2_int) > 0)][which(temp_count_1*temp_count_2 > 0)]

      # r1_id <- c(r1_id, temp_id_r1)
      # r2_id <- c(r2_id, temp_id_r2)

      temp_qa_r1 <- fq_1@quality@quality[which((r1_int*r2_int) > 0)][which(temp_count_1*temp_count_2 > 0)]
      temp_qa_r2 <- fq_2@quality@quality[which((r1_int*r2_int) > 0)][which(temp_count_1*temp_count_2 > 0)]

      temp_q1 <- Biostrings::subseq(temp_qa_r1, delete_1 + 1)
      temp_q2 <- Biostrings::subseq(temp_qa_r2, delete_2 + 1)

      # r1_qa <- c(r1_qa, temp_q1)
      # r2_qa <- c(r2_qa, temp_q2)

      names(temp_r1) <- as.character(temp_id_r1)

      names(temp_r2) <- as.character(temp_id_r2)

      Biostrings::writeXStringSet(temp_r1,
                                  r1_names[i],
                                  format = "fastq",
                                  qualities = temp_q1,
                                  append = TRUE)
      
      Biostrings::writeXStringSet(temp_r2,
                                  r2_names[i],
                                  format = "fastq",
                                  qualities = temp_q2,
                                  append = TRUE)

      cat(paste("     tags: ",
                paste(fwd, rvr, sep = "-"),
                " tag r1: ",
                length(temp_r1),
                " tag r2: ",
                length(temp_r2),
                "\n",
                sep = ""))

    }

  }

  cat("closing connections...\n")
  close(strm_r1)
  close(strm_r2)

}



