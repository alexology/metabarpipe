### TAKEN AS THEY ARE FROM https://github.com/benjjneb/ITS-Workflow/blob/master/ITS_workflow.md
### KUDOS TO benjjneb FOR HIS WONDERFUL WORK
#' @importFrom Biostrings DNAString complement reverse reverseComplement vcountPattern DNAStringSet

prmr_funcn <- function(prmr) {
  prmr_dna <- Biostrings::DNAString(prmr)
  cmplmnt <- toString(Biostrings::complement(prmr_dna))
  rvrse <- toString(Biostrings::reverse(prmr_dna))
  rvrse_cmplmnt <- toString(Biostrings::reverseComplement(prmr_dna))
  return(c(prmr, cmplmnt, rvrse, rvrse_cmplmnt))
}


vcnt_funcn <- function(prmr, filt_seq) {
  vcount_pattern <- Biostrings::vcountPattern(prmr,
                                              Biostrings::DNAStringSet(filt_seq),
                                              fixed = FALSE)
  return(sum(vcount_pattern))
}

# https://stackoverflow.com/questions/12403312/find-the-number-of-spaces-in-a-string
countSpaces <- function(s){sapply(gregexpr(" ", s), function(p) {sum(p>=0)})}
