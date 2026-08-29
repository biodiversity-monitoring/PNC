#' COMBINE Mammal Trait Data
#'
#' @description
#' A processed subset of the COMBINE database containing 40 quantitative
#' mammal traits and ecological attributes. The dataset contains one record
#' for each of 6,234 extant or recently extinct mammal species.
#'
#' @format A data frame with 6,234 rows and 43 variables:
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
#'   \item{adult_mass_g}{
#'     Numeric. Adult body mass, g.
#'   }
#'   \item{adult_brain_mass_g}{
#'     Numeric. Adult brain mass, g.
#'   }
#'   \item{adult_body_length_mm}{
#'     Numeric. Adult body length from the tip of the nose to the anus or
#'     base of the tail, mm.
#'   }
#'   \item{adult_forearm_length_mm}{
#'     Numeric. Adult forearm length from the elbow to the wrist, mm.
#'     This variable applies to bats.
#'   }
#'   \item{max_longevity_d}{
#'     Numeric. Maximum reported lifespan, days.
#'   }
#'   \item{maturity_d}{
#'     Numeric. Age at sexual maturity, days.
#'   }
#'   \item{female_maturity_d}{
#'     Numeric. Age at which a female reaches sexual maturity, days.
#'   }
#'   \item{male_maturity_d}{
#'     Numeric. Age at which a male reaches sexual maturity, days.
#'   }
#'   \item{age_first_reproduction_d}{
#'     Numeric. Age at which a female first gives birth or the young first
#'     attach to the teats, days.
#'   }
#'   \item{gestation_length_d}{
#'     Numeric. Duration of gestation, days.
#'   }
#'   \item{teat_number_n}{
#'     Numeric. Number of teats.
#'   }
#'   \item{litter_size_n}{
#'     Numeric. Number of offspring per litter per female.
#'   }
#'   \item{litters_per_year_n}{
#'     Numeric. Number of litters per female per year.
#'   }
#'   \item{interbirth_interval_d}{
#'     Numeric. Time between successive reproductive events, days.
#'   }
#'   \item{neonate_mass_g}{
#'     Numeric. Body mass at birth, g.
#'   }
#'   \item{weaning_age_d}{
#'     Numeric. Age at which primary nutritional dependence on the mother
#'     ends, days.
#'   }
#'   \item{weaning_mass_g}{
#'     Numeric. Body mass at weaning, g.
#'   }
#'   \item{generation_length_d}{
#'     Numeric. Mean age of the parents of the current cohort, days.
#'   }
#'   \item{dispersal_km}{
#'     Numeric. Distance between the natal location and the location of
#'     reproduction, km.
#'   }
#'   \item{density_n_km2}{
#'     Numeric. Population density, individuals per square kilometre.
#'   }
#'   \item{home_range_km2}{
#'     Numeric. Area used for the routine activities of an individual or
#'     social group, square kilometres.
#'   }
#'   \item{social_group_n}{
#'     Numeric. Number of individuals in a social group.
#'   }
#'   \item{dphy_invertebrate}{
#'     Numeric. Broad diet fraction composed of invertebrates, percent.
#'   }
#'   \item{dphy_vertebrate}{
#'     Numeric. Broad diet fraction composed of vertebrates, percent.
#'   }
#'   \item{dphy_plant}{
#'     Numeric. Broad diet fraction composed of plants or fungi, percent.
#'   }
#'   \item{det_inv}{
#'     Numeric. Detailed diet fraction composed of invertebrates, percent.
#'   }
#'   \item{det_vend}{
#'     Numeric. Detailed diet fraction composed of endothermic vertebrates,
#'     including mammals and birds, percent.
#'   }
#'   \item{det_vect}{
#'     Numeric. Detailed diet fraction composed of ectothermic vertebrates,
#'     including reptiles and amphibians, percent.
#'   }
#'   \item{det_vfish}{
#'     Numeric. Detailed diet fraction composed of fish, percent.
#'   }
#'   \item{det_vunk}{
#'     Numeric. Detailed diet fraction composed of unspecified vertebrates,
#'     percent.
#'   }
#'   \item{det_scav}{
#'     Numeric. Detailed diet fraction composed of scavenged material or
#'     carrion, percent.
#'   }
#'   \item{det_fruit}{
#'     Numeric. Detailed diet fraction composed of fruits or drupes, percent.
#'   }
#'   \item{det_nect}{
#'     Numeric. Detailed diet fraction composed of nectar, pollen, or plant
#'     exudates, percent.
#'   }
#'   \item{det_seed}{
#'     Numeric. Detailed diet fraction composed of seeds, grains, nuts, or
#'     spores, percent.
#'   }
#'   \item{det_plantother}{
#'     Numeric. Detailed diet fraction composed of other plant material,
#'     percent.
#'   }
#'   \item{det_diet_breadth_n}{
#'     Numeric. Number of detailed EltonTraits diet categories that each
#'     contribute at least 20 percent of the diet.
#'   }
#'   \item{upper_elevation_m}{
#'     Numeric. Upper known elevational limit, m.
#'   }
#'   \item{lower_elevation_m}{
#'     Numeric. Lower known elevational limit, m.
#'   }
#'   \item{altitude_breadth_m}{
#'     Numeric. Difference between the upper and lower elevational limits, m.
#'   }
#'   \item{habitat_breadth_n}{
#'     Numeric. Number of distinct suitable level 1 IUCN habitat categories.
#'   }
#' }
#'
#' @details
#' Species names, trait values, and units follow the COMBINE database. This
#' processed dataset contains one record per species. No additional unit
#' conversions, trait transformations, or imputations were applied.
#'
#' @source
#' COMBINE data repository:
#' \doi{10.6084/m9.figshare.13028255.v4}
#'
#' @references
#' Soria, C. D., Pacifici, M., Di Marco, M., Stephen, S. M., and
#' Rondinini, C. (2021). COMBINE: a coalesced mammal database of
#' intrinsic and extrinsic traits. \emph{Ecology}, 102, e03344.
#' \doi{10.1002/ecy.3344}
#'
#' Users of individual traits should also consult and cite the relevant
#' underlying data sources listed in the COMBINE data repository.
#'
#' @examples
#' data(COMBINE)
#'
#' dim(COMBINE)
#' head(COMBINE[, c("species", "adult_mass_g", "maturity_d")])
#'
"COMBINE"
