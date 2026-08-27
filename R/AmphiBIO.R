#' AmphiBIO Amphibian Trait Data
#'
#' @description
#' A processed subset of the AmphiBIO database containing 12 quantitative
#' morphological, life-history, and reproductive traits. The dataset contains
#' one record for each of 6,776 amphibian species.
#'
#' @format A data frame with 6,776 rows and 15 variables:
#' \describe{
#'   \item{species}{
#'     Character. Full binomial species name.
#'   }
#'   \item{genus}{
#'     Character. Genus name.
#'   }
#'   \item{family}{
#'     Character. Taxonomic family.
#'   }
#'   \item{Body_mass_g}{
#'     Numeric. Maximum adult body mass, g.
#'   }
#'   \item{Age_at_maturity_min_y}{
#'     Numeric. Minimum age at maturation or sexual maturity, years.
#'   }
#'   \item{Age_at_maturity_max_y}{
#'     Numeric. Maximum age at maturation or sexual maturity, years.
#'   }
#'   \item{Body_size_mm}{
#'     Numeric. Maximum adult body size, mm.
#'   }
#'   \item{Size_at_maturity_min_mm}{
#'     Numeric. Minimum body size at maturation or sexual maturity, mm.
#'   }
#'   \item{Size_at_maturity_max_mm}{
#'     Numeric. Maximum body size at maturation or sexual maturity, mm.
#'   }
#'   \item{Longevity_max_y}{
#'     Numeric. Maximum lifespan, years.
#'   }
#'   \item{Litter_size_min_n}{
#'     Numeric. Minimum number of offspring or eggs per clutch.
#'   }
#'   \item{Litter_size_max_n}{
#'     Numeric. Maximum number of offspring or eggs per clutch.
#'   }
#'   \item{Reproductive_output_y}{
#'     Numeric. Maximum number of reproductive events per year.
#'   }
#'   \item{Offspring_size_min_mm}{
#'     Numeric. Minimum offspring or egg size, mm.
#'   }
#'   \item{Offspring_size_max_mm}{
#'     Numeric. Maximum offspring or egg size, mm.
#'   }
#' }
#'
#' @details
#' The source AmphiBIO database contains 17 amphibian ecological traits. This
#' processed dataset retains 12 numeric variables selected for use with PNC.
#' These variables include continuous measurements and count variables.
#'
#' Species names and higher taxonomy were retained from AmphiBIO. The source
#' taxonomy follows Amphibian Species of the World version 5.5. Trait values
#' and units were retained as provided in the source database. No additional
#' unit conversions, transformations, or imputations were performed. Missing
#' values are represented by \code{NA}.
#'
#' Maximum body size is generally reported as snout-vent length for Anura and
#' as total length for Gymnophiona and Caudata. For some Caudata records, only
#' snout-vent length was available in the source database.
#'
#' Offspring size generally refers to egg-yolk diameter, although some source
#' studies measured total egg diameter including the surrounding jelly
#' capsule. The measurement convention cannot always be distinguished from the
#' processed data.
#'
#' Trait coverage varies substantially among variables. In this dataset,
#' 57 species have complete values for all 12 selected traits, whereas
#' 796 species have no observed values for any of the selected traits.
#' Missing-data coverage and the dependence between paired minimum and maximum
#' variables should therefore be considered when selecting variables for
#' multivariate analyses.
#'
#' @source
#' AmphiBIO data repository:
#' \url{https://doi.org/10.6084/m9.figshare.4644424}
#'
#' @references
#' Oliveira, B. F., Sao-Pedro, V. A., Santos-Barrera, G., Penone, C., and
#' Costa, G. C. (2017). AmphiBIO, a global database for amphibian ecological
#' traits. \emph{Scientific Data}, 4, 170123.
#' \doi{10.1038/sdata.2017.123}
#'
#' @examples
#' data(AmphiBIO)
#'
#' dim(AmphiBIO)
#' head(AmphiBIO[, c(
#'   "species",
#'   "Body_size_mm",
#'   "Reproductive_output_y"
#' )])
#'
"AmphiBIO"
