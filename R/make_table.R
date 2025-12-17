#' Make a fastq table
#'
#' @description
#' This function make a sample table from e ENA table.
#'
#' @param x A table downloaded from ENA.
#'
#' @importFrom dplyr select slice pull
#' @importFrom httr content GET add_headers
#' @importFrom xml2 as_list
#' @importFrom tibble as_tibble

make_table <- function(x){

  # set the data.frame to store the results
  final_df <-  data.frame(fastq = character(),
                          original = character())

  for(i in 1:nrow(x)){
    # set the first part of the request
    a <- "https://www.ebi.ac.uk/ena/browser/api/xml/"

    # get the accession name
    b <- x %>%
      dplyr::select(1) %>%
      dplyr::slice(i) %>%
      dplyr::pull()

    # set the last part of the request
    d <- "?download=false&gzip=false&includeLinks=false"

    # make the request
    res <- httr::GET(paste(a, b, d, sep = ""),
                     httr::add_headers('Accept: application/json'))

    # transform the request to list
    z <- res %>%
      httr::content(., encoding = "UTF-8") %>%
      xml2::as_list()

    # extract data to populate the data.frame
    final_df[i, 1] <- b
    final_df[i, 2] <- z$RUN_SET$RUN$EXPERIMENT_REF$IDENTIFIERS$SUBMITTER_ID[[1]]

    # print the progression
    message(paste(b, ": ", i, " of ", nrow(x),  sep = "") )
  }

  # output the results
  tibble::as_tibble(final_df)
}
