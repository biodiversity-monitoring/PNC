#' Plant Functional Traits from the TRY Database
#'
#' A species-level compilation of 20 quantitative plant functional traits
#' obtained from the TRY Plant Trait Database. Taxonomic names were standardized
#' using the \pkg{U.Taxonstand} package. Each row represents one species, while
#' the \code{genus} and \code{family} columns provide its higher-taxonomic
#' classification.
#'
#' @format A data frame with 58,975 rows and 23 variables:
#' \describe{
#'   \item{species}{Character. Standardized species name.}
#'   \item{genus}{Character. Standardized genus name.}
#'   \item{family}{Character. Standardized family name.}
#'   \item{DispersalUnitLength}{
#'     Numeric. Dispersal unit length, mm. (TraitID: 237)
#'   }
#'   \item{LA}{
#'     Numeric. Leaf area, or leaflet area for compound leaves, mm^2.
#'     Inclusion of the petiole may vary among source records. (TraitID: 3113)
#'   }
#'   \item{LDMC}{
#'     Numeric. Leaf dry mass per leaf fresh mass, or leaf dry matter content,
#'     g g^-1. (TraitID: 47)
#'   }
#'   \item{LeafC}{
#'     Numeric. Leaf carbon content per unit leaf dry mass, mg g^-1.
#'     (TraitID: 13)
#'   }
#'   \item{LeafN}{
#'     Numeric. Leaf nitrogen content per unit leaf dry mass, mg g^-1.
#'     (TraitID: 14)
#'   }
#'   \item{LeafNPratio}{
#'     Numeric. Leaf nitrogen-to-phosphorus mass ratio, g g^-1
#'     (dimensionless). (TraitID: 56)
#'   }
#'   \item{LeafNperArea}{
#'     Numeric. Leaf nitrogen content per unit leaf area, g m^-2.
#'     (TraitID: 50)
#'   }
#'   \item{LeafP}{
#'     Numeric. Leaf phosphorus content per unit leaf dry mass, mg g^-1.
#'     (TraitID: 15)
#'   }
#'   \item{Leafdelta15N}{
#'     Numeric. Leaf nitrogen isotope signature (delta 15N), per mille
#'     (‰). (TraitID: 78)
#'   }
#'   \item{Leaffreshmass}{
#'     Numeric. Leaf fresh mass, g. (TraitID: 163)
#'   }
#'   \item{LMA}{
#'     Numeric. Leaf mass per area, calculated as the reciprocal of
#'     \code{SLA}, mg mm^-2.
#'   }
#'   \item{PlantHeight}{
#'     Numeric. Vegetative plant height, m. (TraitID: 3106)
#'   }
#'   \item{RootingDepth}{
#'     Numeric. Rooting depth, m. (TraitID: 6)
#'   }
#'   \item{SeedLength}{
#'     Numeric. Seed length, mm. (TraitID: 27)
#'   }
#'   \item{SeedMass}{
#'     Numeric. Seed dry mass, mg. (TraitID: 26)
#'   }
#'   \item{SeedNumber}{
#'     Numeric. Number of seeds per reproductive unit. (TraitID: 138)
#'   }
#'   \item{SLA}{
#'     Numeric. Leaf area per unit leaf dry mass, or specific leaf area,
#'     with the petiole excluded, mm^2 mg^-1. (TraitID: 3115)
#'   }
#'   \item{SSD}{
#'     Numeric. Stem specific density, defined as stem dry mass per unit
#'     fresh volume, or wood density, g cm^-3. (TraitID: 4)
#'   }
#'   \item{StemConduitDensity}{
#'     Numeric. Density of stem conduits, including vessels and tracheids,
#'     mm^-2. (TraitID: 169)
#'   }
#'   \item{WoodVesselLength}{
#'     Numeric. Length of stem conduit elements, including vessel elements
#'     and tracheids, µm. (TraitID: 282)
#'   }
#' }
#'
#' @details
#' The dataset contains 19 traits obtained directly from TRY and one derived
#' variable, \code{LMA}. The latter was calculated as the reciprocal of
#' \code{SLA} and therefore has no separate TRY TraitID. Because \code{SLA}
#' is expressed as mm^2 mg^-1, \code{LMA} is expressed as mg mm^-2 in this
#' dataset.
#'
#' The traits describe major dimensions of leaf structure and chemistry,
#' plant size and rooting strategy, seed and dispersal strategy, and stem
#' and wood anatomy. Missing trait values are represented by \code{NA}.
#'
#' @source
#' TRY Plant Trait Database (\url{https://www.try-db.org/})
#'
#' @references
#' Kattge, J., Bönisch, G., Díaz, S., et al. (2020). TRY plant trait
#' database: enhanced coverage and open access. Global Change Biology,
#' 26, 119–188. \doi{10.1111/gcb.14904}
#'
#' Zhang, J., and Qian, H. (2023). U.Taxonstand: An R package for
#' standardizing scientific names of plants and animals. Plant Diversity,
#' 45, 1–5. \doi{10.1016/j.pld.2022.09.001}
#'
#' @examples
#' data("TRY")
#' head(TRY)
#'
#' @keywords datasets
"TRY"
