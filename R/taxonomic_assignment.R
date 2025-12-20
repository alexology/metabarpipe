#' taxonomic assignment
#' 
#' @details
#' This function build a table using boldigger dta
#' 
#' @param project_path Path where searching for boldigger results
#' @param asv_folder Path to the folder with ASV results. 
#' @param species_only Keep records with species name only.
#' @param remove_bad_names Remove badly formatted species names.
#' @param filter_taxlev Select a taxonomic level to subset (e.g. Phylum)
#' @param filter_taxa Select one taxon or more taxa to filter at the selected taxonomic level.
#'
#' @export
#' 
#' 
#'
#' @importFrom readxl read_excel
#' @importFrom dplyr all_of select filter inner_join pull
#' @importFrom tidyr unite separate
#' @importFrom Biostrings readDNAStringSet
#' @importFrom writexl write_xlsx
#' @importFrom tibble as_tibble rownames_to_column
#' @importFrom rlang .data



taxonomic_assignment <- function(project_path = NULL,
                                 asv_folder = NULL,
                                 species_only = FALSE,
                                 remove_bad_names = TRUE,
                                 filter_taxlev = NULL,
                                 filter_taxa = NULL){
  
  if(is.null(project_path)){
    project_path <- file.path(getwd())
  }
  
  # colnames lower case for aaplotype 
  if(! is.null(filter_taxlev)){
    filter_taxlev <- tolower(filter_taxlev)
  }
  
  # list asv file
  asv_file <- list.files(file.path(project_path, asv_folder), full.names = TRUE)
  
  # read the asv file
  asv_file <- readxl::read_excel(asv_file)
  
  # get boldigger results
  boldigger_results <- readxl::read_excel(file.path(project_path,
                                                    "13_taxonomic_assignment",
                                                    "boldigger3_data",
                                                    "boldigger_identification_result.xlsx"))
  
  # lower case column for aaplotype
  colnames(boldigger_results) <- tolower(colnames(boldigger_results))
  
  # boldigger fasta, needed for merging
  boldigger_results <- Biostrings::readDNAStringSet(file.path(project_path,
                                                            "13_taxonomic_assignment",
                                                            "boldigger.fasta")) %>%
    as.data.frame() %>%
    dplyr::rename("ASV" = 1) %>%
    tibble::rownames_to_column("id") %>%
    dplyr::inner_join(boldigger_results, by = "id") %>%
    dplyr::select(-c("id") ) %>%
    tibble::as_tibble()
  
  if(isTRUE(remove_bad_names)){
    boldigger_results <- boldigger_results %>%
      dplyr::filter(!grepl("cf.", .data$species, fixed = TRUE)) %>%
      dplyr::filter(!grepl("sp.", .data$species, fixed = TRUE)) %>%
      dplyr::filter(!grepl("sp.", .data$species, fixed = TRUE)) %>%
      dplyr::filter(!grepl("sp.", .data$species, fixed = TRUE)) %>%
      dplyr::filter(!grepl("/", .data$species, fixed = TRUE)) %>%
      dplyr::filter(!grepl("no-match", .data$species, fixed = TRUE)) %>%
      dplyr::filter(!grepl(".*[0-9].*", .data$species)) 
    
  }
  
  if(isTRUE(species_only)){
    boldigger_results <- boldigger_results %>%
      dplyr::filter(countSpaces(.data$species) == 1) %>%
      dplyr::filter(.data$species != "") %>%
      dplyr::filter(!is.na(.data$species)) %>%
      dplyr::filter(!is.na(.data$genus))
  } 

  if(is.null(filter_taxlev) & !is.null(filter_taxa)){
    if(!is.null(filter_taxa)){
      stop("filter_taxlev cannot be NULL when filter_taxa is not NULL")
    }
  }
    
  
  if(!is.null(filter_taxlev)){
    if(is.null(filter_taxa)){
      stop("filter_taxa cannot be NULL when filter_taxlev is not NULL")
    }
    
    # subset a taxon as requested by the user
    #boldigger_results <- boldigger_results[dplyr::pull(boldigger_results[, filter_taxlev]) %in% filter_taxa,]
    boldigger_results <- boldigger_results %>%
      dplyr::filter(.data[[filter_taxlev]] %in% filter_taxa)
  }
  
  boldigger_results <- boldigger_results %>%
    dplyr::select(dplyr::all_of("ASV"):dplyr::all_of("species")) %>%
    dplyr::inner_join(asv_file, by = "ASV")
  
  # set the name of the output to keep trace of the calculations
  file_name <- c("")
 
  if(isTRUE(remove_bad_names)){
    file_name <- paste(file_name, "bad_removed", sep = "_")
  }
   
  if(isTRUE(species_only)){
    file_name <- paste(file_name, "species_only", sep = "_")
  }
  
  if(!is.null(filter_taxa)){
    file_name <- paste(file_name, "taxa_filtered", sep = "_")
  }
  
  file_name <- paste("taxonomic_assigned", file_name, ".xlsx", sep = "")
  
  
  
  writexl::write_xlsx(boldigger_results, file.path(project_path, "13_taxonomic_assignment", file_name))
  boldigger_results  
}