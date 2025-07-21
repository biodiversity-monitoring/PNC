#' Extract Plant Traits from Trait Database
#'
#' This function extracts plant trait data from the TRY database or similar datasets
#' for a specified list of taxa at different taxonomic ranks (species, genus, or family).
#' For numeric traits at genus and family levels, it calculates mean values across
#' all available records.
#'
#' @param sp.list A character vector containing the names of taxa to extract traits for.
#'   The names should match the taxonomic rank specified in the 'rank' parameter.
#' @param dataset A data frame containing trait data. Default is TRY database.
#'   Must contain columns named "species", "genus", and "family" for taxonomic information.
#' @param rank A character string specifying the taxonomic rank to match against.
#'   Must be one of "species", "genus", or "family". Default is "species".
#' @param traits A character vector specifying which traits to extract. If NULL (default),
#'   all available traits in the dataset will be extracted. Available traits are all
#'   columns except "species", "genus", and "family".
#'
#' @return A data frame with taxa names as row names and trait names as column names.
#'   For species-level extraction, returns the first occurrence of each species.
#'   For genus/family-level extraction, returns mean values for numeric traits
#'   and the first occurrence for non-numeric traits. Missing values are represented as NA.
#'
#' @details
#' The function performs the following operations:
#' \itemize{
#'   \item Validates input parameters
#'   \item Identifies available traits in the dataset
#'   \item Matches input taxa with dataset entries
#'   \item Reports missing taxa
#'   \item Extracts trait data based on the specified taxonomic rank
#'   \item For numeric traits at genus/family level, calculates mean values
#'   \item For non-numeric traits, uses the first available value
#'   \item Handles NaN values by converting them to NA
#' }
#'
#' @examples
#' # Load the dataset
#' data(TRY)
#'
#' # Extract all traits for species
#' species_list <- c("Acaena novae-zelandiae", "Adiantum capillus-veneris", "Zuelania guidonia")
#' extract_traits(species_list, TRY, rank = "species")
#'
#' # Extract specific traits for species
#' extract_traits(species_list, TRY, rank = "species",
#'                traits = c("LA", "LMA", "LeafN", "PlantHeight", "SeedMass", "SSD"))
#'
#' # Extract specific traits at genus level
#' genus_list <- c("Acaena", "Adiantum")
#' extract_traits(genus_list, TRY, rank = "genus",
#'                traits = c("LDMC", "PlantHeight", "SeedMass"))
#'
#' @export
extract_traits <- function(sp.list, dataset, rank = "species", traits = NULL) {
  if (!rank %in% c("species", "genus", "family")) {
    stop("rank must be one of 'species', 'genus', or 'family'")
  }
  if (is.null(traits)) {
    trait_cols <- base::setdiff(base::colnames(dataset), c("species", "genus", "family"))
  } else {
    available_traits <- base::setdiff(base::colnames(dataset), c("species", "genus", "family"))
    invalid_traits <- base::setdiff(traits, available_traits)
    if (base::length(invalid_traits) > 0) {
      base::cat("Warning: The following traits are not available in the dataset:\n")
      base::print(invalid_traits)
    }
    trait_cols <- base::intersect(traits, available_traits)
    if (base::length(trait_cols) == 0) {
      stop("No valid traits specified")
    }
  }
  result <- base::data.frame(
    base::matrix(NA, nrow = base::length(sp.list), ncol = base::length(trait_cols))
  )
  base::rownames(result) <- sp.list
  base::colnames(result) <- trait_cols
  trait_taxa <- dataset[[rank]]
  matched_taxa <- base::intersect(sp.list, trait_taxa)
  missing_taxa <- base::setdiff(sp.list, trait_taxa)
  base::cat("The missing", rank, "in the dataset: (", base::length(missing_taxa), "):\n")
  base::print(missing_taxa)
  if (base::length(matched_taxa) > 0) {
    for (taxon in matched_taxa) {
      taxon_rows <- base::which(dataset[[rank]] == taxon)
      if (base::length(taxon_rows) > 0) {
        if (rank == "species") {
          result[taxon, ] <- dataset[taxon_rows[1], trait_cols]
        } else {
          taxon_data <- dataset[taxon_rows, ]
          for (col in trait_cols) {
            if (base::is.numeric(taxon_data[[col]])) {
              col_mean <- base::mean(taxon_data[[col]], na.rm = TRUE)
              result[taxon, col] <- if (base::is.nan(col_mean)) NA else col_mean
            } else {
              result[taxon, col] <- taxon_data[[col]][1]
            }
          }
        }
      }
    }
  }
  result[] <- base::lapply(result, function(x) {
    x[base::is.nan(x)] <- NA
    x
  })
  return(result)
}
