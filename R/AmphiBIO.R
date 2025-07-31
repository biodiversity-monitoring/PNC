#' AmphiBIO: Global Amphibian Ecological Traits Database
#'
#' @description
#' A comprehensive global database of ecological traits for amphibian species,
#' compiled to provide insights into the life history and ecological characteristics
#' of amphibians worldwide.
#'
#' @format A data frame with multiple variables:
#' \describe{
#'   \item{species}{Scientific name of the amphibian species}
#'   \item{genus}{Taxonomic genus of the species}
#'   \item{family}{Taxonomic family of the species}
#'   \item{Body_mass_g}{Maximum adult body mass.}
#'   \item{Age_at_maturity_min_y}{Minimum age at maturation or sexual maturity.}
#'   \item{Age_at_maturity_max_y}{Maximum age at maturation or sexual maturity.}
#'   \item{Body_size_mm}{Maximum adult body size. In Anura, body size is reported as snout to vent length. In Gymnophiona and Caudata, body size is reported as total length.}
#'   \item{Size_at_maturity_min_mm}{Minimum size at maturation or sexual maturity.}
#'   \item{Size_at_maturity_max_mm}{Maximum size at maturation or sexual maturity.}
#'   \item{Longevity_max_y}{Maximum life span.}
#'   \item{Litter_size_min_n}{Minimum no. of offspring or eggs per clutch.}
#'   \item{Litter_size_max_n}{Maximum no. of offspring or eggs per clutch.}
#'   \item{Reproductive_output_y}{Maximum no. reproduction events per year.}
#'   \item{Offspring_size_min_mm}{Minimum offspring or egg size.}
#'   \item{Offspring_size_max_mm}{Maximum offspring or egg size.}
#' }
#'
#' @references
#' Oliveira, B. F., São-Pedro, V. A., Santos-Barrera, G., Penone, C., & Costa, G. C. (2017). AmphiBIO, a global database for amphibian ecological traits. Scientific data, 4(1), 1-7.
#' \url{https://doi.org/10.1038/sdata.2017.123}
#'
#' @examples
#' # Load the dataset
#' data(AmphiBIO)
#' head(AmphiBIO)
"AmphiBIO"
