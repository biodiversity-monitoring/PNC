#' ReptTraits: A Comprehensive Dataset of Ecological Traits in Reptiles
#'
#' A comprehensive dataset containing ecological and morphological characteristics of reptiles.
#' The dataset provides detailed information about reptile species, including elevation,
#' seasonal precipitation, body mass, and reproductive features.
#'
#' @format A data frame with the following columns:
#' \describe{
#'   \item{species}{Scientific species name}
#'   \item{genus}{Genus name}
#'   \item{family}{Family name}
#'   \item{Minimal_elevation}{Minimum elevation where the species was observed (meters above sea level)}
#'   \item{Maximum_elevation}{Maximum elevation where the species was observed (meters above sea level)}
#'   \item{Seasonality_Precipitation}{Seasonal precipitation information}
#'   \item{Maximum_body_mass}{Maximum body mass of the species (grams)}
#'   \item{Maximum_length}{Maximum length ("SVL", mm)/straight carapace length for turtles ("SCL", mm)}
#'   \item{Mean_number_of_offspring}{Mean number of offspring or eggs per clutch}
#'   \item{Smallest_clutch_size}{Minimum clutch/litter size}
#'   \item{Largest_clutch_size}{Maximum clutch/litter size}
#' }
#'
#' @references
#' Oskyrko, O., Mi, C., Meiri, S., & Du, W. (2024). ReptTraits: a comprehensive dataset of ecological traits in reptiles. Scientific Data, 11(1), 243.
#' \url{https://doi.org/10.1038/s41597-024-03079-5}
#'
#' @examples
#' data(ReptTraits)
#' head(ReptTraits)
"ReptTraits"
