#' Bird Morphological Traits from AVONET
#'
#' A species-level compilation of 11 quantitative morphological traits
#' obtained from AVONET. The dataset also includes genus and family
#' information based on BirdLife International taxonomy. Each row
#' represents one bird species.
#'
#' @format A data frame with 11,009 rows and 14 variables:
#' \describe{
#'   \item{species}{
#'     Character. Scientific species name.
#'   }
#'   \item{genus}{
#'     Character. Genus name.
#'   }
#'   \item{family}{
#'     Character. Family name based on BirdLife International taxonomy.
#'   }
#'   \item{Beak.Length_Culmen}{
#'     Numeric. Beak length measured along the culmen from the beak tip
#'     to the base of the skull, mm.
#'   }
#'   \item{Beak.Length_Nares}{
#'     Numeric. Beak length from the anterior edge of the nares to the
#'     beak tip, mm.
#'   }
#'   \item{Beak.Width}{
#'     Numeric. Beak width at the anterior edge of the nares, mm.
#'   }
#'   \item{Beak.Depth}{
#'     Numeric. Beak depth at the anterior edge of the nares, mm.
#'   }
#'   \item{Tarsus.Length}{
#'     Numeric. Tarsus length from the posterior notch between the tibia
#'     and tarsus to the end of the last scale of the acrotarsium, mm.
#'   }
#'   \item{Wing.Length}{
#'     Numeric. Length from the carpal joint to the tip of the longest
#'     primary feather on the unflattened wing, mm.
#'   }
#'   \item{Kipps.Distance}{
#'     Numeric. Distance between the tip of the first outermost secondary
#'     feather and the tip of the longest primary feather, mm. This value
#'     was measured directly or calculated as \code{Wing.Length -
#'     Secondary1}.
#'   }
#'   \item{Secondary1}{
#'     Numeric. Length from the carpal joint to the tip of the first
#'     outermost secondary feather, mm.
#'   }
#'   \item{Hand-Wing.Index}{
#'     Numeric. Hand-wing index, calculated as
#'     \code{100 * Kipps.Distance / Wing.Length}; dimensionless.
#'   }
#'   \item{Tail.Length}{
#'     Numeric. Distance from the tip of the longest rectrix to the point
#'     where the two central rectrices emerge from the skin, mm.
#'   }
#'   \item{Mass}{
#'     Numeric. Species-level mean body mass across available individuals,
#'     including both sexes, g.
#'   }
#' }
#'
#' @details
#' Values are species-level summaries provided by AVONET. The four beak
#' measurements describe feeding morphology, whereas measurements of the
#' tarsus, wing, secondary feathers, and tail describe locomotor morphology.
#' The hand-wing index measures wing elongation and is commonly used as an
#' indicator of flight efficiency and dispersal ability. Body mass represents
#' overall body size.
#'
#' This PNC data object contains the 11 morphological traits and associated
#' taxonomic information. The ecological and geographical variables included
#' in the complete AVONET release are not included. Missing values are
#' represented by \code{NA}.
#'
#' @source
#' AVONET data repository
#' (\url{https://figshare.com/s/b990722d72a26b5bfead})
#'
#' @references
#' Tobias, J. A., Sheard, C., Pigot, A. L., et al. (2022).
#' AVONET: morphological, ecological and geographical data for all birds.
#' Ecology Letters, 25, 581–597. \doi{10.1111/ele.13898}
#'
#' @examples
#' data("AVONET")
#' head(AVONET)
#'
#' @keywords datasets
"AVONET"
