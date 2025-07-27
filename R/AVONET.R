#' AVONET Bird Morphological Dataset
#'
#' @description
#' Comprehensive morphological dataset for bird species, including taxonomic information
#' from BirdLife International and detailed morphological measurements.
#'
#' @format A data frame with 11,009 rows and 14 columns, where each row represents a bird species:
#' \describe{
#'   \item{species}{Species scientific name}
#'   \item{genus}{Genus name}
#'   \item{family}{Family name, according to BirdLife International taxonomy}
#'   \item{Beak.Length_Culmen}{Length from beak tip to skull base, in millimeters}
#'   \item{Beak.Length_Nares}{Length from nostril anterior edge to beak tip, in millimeters}
#'   \item{Beak.Width}{Beak width at the anterior edge of nostrils, in millimeters}
#'   \item{Beak.Depth}{Beak depth at the anterior edge of nostrils, in millimeters}
#'   \item{Tarsus.Length}{Tarsus length from posterior notch between tibia and tarsus to the last scale end, in millimeters}
#'   \item{Wing.Length}{Length from carpal joint to longest primary feather tip, in millimeters}
#'   \item{Kipps.Distance}{Length from first secondary feather tip to longest primary feather tip, in millimeters}
#'   \item{Secondary1}{Length from carpal joint to first secondary feather tip, in millimeters}
#'   \item{Hand-Wing.Index}{100*DK/Lw, where DK is Kipp's distance and Lw is wing length}
#'   \item{Tail.Length}{Distance from longest rectrix tip to point where central rectrices protrude from skin, in millimeters}
#'   \item{Mass}{Species average body mass, including both male and female, in grams}
#' }
#'
#' @details
#' This dataset provides comprehensive morphological measurements of birds,
#' including beak, wing, tarsus, and body weight indicators.
#' Data originates from a comprehensive study of bird morphological,
#' ecological, and geographical characteristics.
#'
#' @note
#' - Taxonomic information based on BirdLife International
#' - Measurements represent species averages
#' - Hand-Wing Index reflects flight capability and ecological adaptation
#'
#' @references
#' Tobias, J.A., Sheard, C., Pigot, A.L., et al. (2022)
#' AVONET: morphological, ecological and geographical data for all birds.
#' Ecology Letters, 25, 581–597.
#' \url{https://doi.org/10.1111/ele.13898}
#'
#' @examples
#' data(AVONET)
#' head(AVONET)
#'
"AVONET"
