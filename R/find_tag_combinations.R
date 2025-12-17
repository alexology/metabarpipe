#' Find tag combinations in multiplexed fastq
#'
#' @description
#' This function tries to infer tags from the multiplexed fastq.
#' 
#' @param r1 Path to the multiplexed r1 file.
#' @param r2 Path to the multiplexed r2 file.
#' @param primer_fwd Forward primers.
#' @param primer_rvr Reverse primers.
#' @param project_path The path to the project folder.
#' 
#' @details
#' This function starts by filtering the multiplexed fastq by an expected
#' error of 3 and a sequence length of 300. The cut the first 6 bases and
#' count the frequency of the strings in the overall files. A summary file
#' and the inferred indices are written to disk. 
#' 
#' @export
#' 
#' @importFrom Biostrings subseq vmatchPattern
#' @importFrom dada2 filterAndTrim
#' @importFrom dplyr arrange as_tibble bind_rows filter left_join rename select pull
#' @importFrom tidyr separate
#' @importFrom ShortRead FastqStreamer yield
#' @importFrom writexl write_xlsx


find_tag_combinations <- function(r1, r2, primer_fwd, primer_rvr, project_path){
  
  # trim by length and quality with data2 
  # get r1 name and extension
  r1_splitted_name <- basename(r1) %>%
    strsplit(., "\\.") %>%
    unlist()
  
  # modify the first element for saving the results
  r1_new_name <- paste(r1_splitted_name[1],
                       "_trimmed.",
                       r1_splitted_name[2],
                       sep ="")
  
  # get the full length path
  r1_new_name <- gsub(basename(r1), r1_new_name, r1)
  
  # get r1 name and extension
  r2_splitted_name <- basename(r2) %>%
    strsplit(., "\\.") %>%
    unlist()
  
  # modify the first element for saving the results
  r2_new_name <- paste(r2_splitted_name[1],
                       "_trimmed.",
                       r2_splitted_name[2],
                       sep ="")
  
  # get the full length path
  r2_new_name <- gsub(basename(r2), r2_new_name, r2)
  
  
  dada2::filterAndTrim(fwd = r1,
                       filt = r1_new_name,
                       rev = r2,
                       filt.rev = r2_new_name,
                       minLen = 200,
                       maxN = 0,
                       maxEE = 3)
  
  # set the function to import chunks of fastq in R
  strm_r1 <- ShortRead::FastqStreamer(r1_new_name)
  strm_r2 <- ShortRead::FastqStreamer(r2_new_name)
  
  # set the data.frame where store the results
  res <- data.frame(tag_1 = character(),
                    tag_2 = character(),
                    combo = character())
  
  indices <- data.frame(barcode = character(),
                        rm = numeric())
  
  repeat {
    
    # load r1 and r2
    fq_1 <- ShortRead::yield(strm_r1)
    
    if(length(fq_1) == 0) {break()}
    
    ## process chunk
    fq_2 <- ShortRead::yield(strm_r2)
    
    # find the position of the forward primer
    r1_int_f <- Biostrings::vmatchPattern(primer_fwd,
                                          fq_1@sread,
                                          fixed = FALSE,
                                          max.mismatch = 0,
    )
    
    # transform the results to data.frame, doing this way rows with no
    # primer match are excluded
    r1_int_f_df <- as.data.frame(r1_int_f)
    
    # include sequences without match
    r1_int_f_df <- dplyr::left_join(data.frame(group = 1:length(fq_1@sread)),
                                    r1_int_f_df,
                                    by = "group")
    
    # take the first match if for a sequences more primer matches are found
    # this assure that the selected match is the nearest to the beginning
    # of the sequence
    r1_int_f_df <- r1_int_f_df[! duplicated(r1_int_f_df$group),] %>%
      dplyr::select(1, 3, 4)
    
    # do the same for the reverse primer if the P5-P7 method is used
    if(!is.null(primer_rvr)){
      # find the position of the reverse primer
      r1_int_r <- Biostrings::vmatchPattern(primer_rvr,
                                            fq_1@sread,
                                            fixed = FALSE,
                                            max.mismatch = 0)
      
      # transform the results to data.frame, doing this way rows with no
      # primer match are excluded
      r1_int_r_df <- as.data.frame(r1_int_r)
      
      r1_int_r_df <- dplyr::left_join(data.frame(group = 1:length(fq_1@sread)),
                                      r1_int_r_df,
                                      by = "group")
      
      r1_int_r_df <- r1_int_r_df[! duplicated(r1_int_r_df$group),] %>%
        dplyr::select(1, 3, 4)
      
      r1_int_f_df$start <- suppressWarnings(apply(data.frame(r1_int_f_df$start, r1_int_r_df$start),
                                                  1,
                                                  function(x) min(x, na.rm = TRUE))) 
      
      r1_int_f_df$end <- suppressWarnings(apply(data.frame(r1_int_f_df$end, r1_int_r_df$end),
                                                1,
                                                function(x) min(x, na.rm = TRUE))) 
      
    }
    
    # r1_int_f_df$start[is.na(r1_int_f_df$start)] <- 0
    r1_int_f_df$end[is.infinite(r1_int_f_df$end)] <- NA
    r1_int_f_df$start[is.infinite(r1_int_f_df$start)] <- NA
    
    r1_index <- r1_int_f_df$end
    r1_index[r1_index < 0] <- NA
    
    r1_tags <- Biostrings::subseq(fq_1@sread, 1,  6) %>%
      as.character()
    
    r1_tags[is.na(r1_index)] <- NA
    
    r1_indices <- data.frame(barcode = r1_tags, 
                             rm = r1_int_f_df$start - 1)
    
    
    
    # find the position of the forward primer
    r2_int_f <- Biostrings::vmatchPattern(primer_fwd,
                                          fq_2@sread,
                                          fixed = FALSE,
                                          max.mismatch = 0
    )
    
    # transform the results to data.frame, doing this way rows with no
    # primer match are excluded
    r2_int_f_df <- as.data.frame(r2_int_f)
    
    
    r2_int_f_df <- dplyr::left_join(data.frame(group = 1:length(fq_2@sread)),
                                    r2_int_f_df,
                                    by = "group")
    
    r2_int_f_df <- r2_int_f_df[! duplicated(r2_int_f_df$group),] %>%
      dplyr::select(1, 3, 4)
    
    
    
    if(!is.null(primer_rvr)){
      # find the position of the reverse primer
      r2_int_r <- Biostrings::vmatchPattern(primer_rvr,
                                            fq_2@sread,
                                            fixed = FALSE,
                                            max.mismatch = 0)
      
      # transform the results to data.frame, doing this way rows with no
      # primer match are excluded
      r2_int_r_df <- as.data.frame(r2_int_r)
      
      r2_int_r_df <- dplyr::left_join(data.frame(group = 1:length(fq_2@sread)),
                                      r2_int_r_df,
                                      by = "group")
      
      r2_int_r_df <- r2_int_r_df[! duplicated(r2_int_r_df$group),] %>%
        dplyr::select(1, 3, 4)
      
      r2_int_f_df$start <- suppressWarnings(apply(data.frame(r2_int_f_df$start, r2_int_r_df$start),
                                                  1,
                                                  function(x) min(x, na.rm = TRUE))) 
      
      r2_int_f_df$end <- suppressWarnings(apply(data.frame(r2_int_f_df$end, r2_int_r_df$end),
                                                1,
                                                function(x) min(x, na.rm = TRUE))) 
    }
    
    r2_int_f_df$end[is.infinite(r2_int_f_df$end)] <- NA
    r2_int_f_df$start[is.infinite(r2_int_f_df$start)] <- NA
    
    r2_index <- r2_int_f_df$end
    r2_index[r2_index < 0] <- NA
    
    r2_tags <- Biostrings::subseq(fq_2@sread, 1,  6) %>%
      as.character()
    
    r2_tags[is.na(r2_index)] <- NA
    
    r2_indices <- data.frame(barcode = r2_tags, 
                             rm = r2_int_f_df$start - 1)
    
    # get ith results
    res_temp <- data.frame(tag_1 = r1_tags,
                           tag_2 = r2_tags,
                           combo = paste(r1_tags, r2_tags, sep = "_"))
    
    # remove rowswith empty strings
    res_temp <- res_temp %>%
      dplyr::filter(!is.na(tag_1)) %>%
      dplyr::filter(!is.na(tag_2))
    
    # bind the results
    res <- rbind(res, res_temp)
    
    # process indices
    indices <- r1_indices %>%
      dplyr::bind_rows(r2_indices) %>%
      na.omit() %>%
      unique() %>%
      dplyr::filter(rm > 0)
    
    # cat
    cat("#")
  }
  
  cat("generating xlsx files\n")
  # arrange the results in the correct format and print them to screen  
  res %>%
    dplyr::pull("combo") %>%
    table() %>%
    as.data.frame() %>%
    dplyr::as_tibble() %>%
    dplyr::arrange(dplyr::desc(Freq)) %>%
    dplyr::rename(combinations = 1) %>%
    tidyr::separate(combinations, c("index_1", "index_2")) %>%
    writexl::write_xlsx(file.path(project_path,
                                  "0_raw_files",
                                  "summary_tags.xlsx"))
  
  cat("summary written to disk\n")
  
  indices %>%
    writexl::write_xlsx(file.path(project_path,
                                  "0_raw_files",
                                  "tag_sequences.xlsx"))
    
  cat("summary written to disk\n")
}
