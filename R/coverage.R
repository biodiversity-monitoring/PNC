#' Summarize Trait Coverage
#'
#' `coverage()` summarizes missingness in a trait dataset for individual
#' traits, complete cases across all supplied traits, and all trait values
#' combined.
#'
#' @param data A data frame or matrix with observations, typically species or
#'   other taxa, in rows and traits in columns. All supplied columns are
#'   treated as traits. Species or taxon identifiers should therefore be stored
#'   as row names or removed before calling the function.
#'
#' @return A data frame with one row for each supplied trait and an additional
#'   `"All"` row. It contains the following columns:
#'   \describe{
#'     \item{Trait}{Trait name. The `"All"` row summarizes complete cases
#'       across all supplied traits.}
#'     \item{Available_count}{For an individual trait, the number of
#'       observations with a non-missing value. For the `"All"` row, the
#'       number of observations with complete data for all traits.}
#'     \item{Missing_count}{For an individual trait, the number of
#'       observations with a missing value. For the `"All"` row, the number
#'       of observations missing at least one trait.}
#'     \item{Trait_coverage_rate}{The percentage of available observations,
#'       returned as a character string containing a percent sign.}
#'   }
#'
#'   Overall coverage across all cells is displayed as a message but is not
#'   included in the returned data frame.
#'
#' @details
#' Coverage is summarized in three ways:
#' \itemize{
#'   \item Trait-specific coverage is the percentage of observations with a
#'     non-missing value for each trait.
#'   \item Coverage of complete cases is the percentage of observations with
#'     non-missing values for every supplied trait and is reported in the
#'     `"All"` row.
#'   \item Overall coverage is the percentage of all cells in the supplied
#'     dataset that contain non-missing values and is displayed when the
#'     function runs.
#' }
#'
#' Missingness is defined using `is.na()`. The function summarizes data
#' availability but does not evaluate whether non-missing values are
#' biologically valid for a particular trait.
#'
#' @examples
#' trait_data <- data.frame(
#'   PlantHeight = c(1.2, 1.5, NA, 2.1, 1.8),
#'   LDMC = c(0.5, NA, 0.8, 1.2, 0.9),
#'   LA = c(15.2, 18.5, 12.3, NA, 16.7)
#' )
#'
#' coverage(trait_data)
#'
#' @export
coverage <- function(data) {

  # ---------------------------------------------------------------------------
  # 1. Check input
  # ---------------------------------------------------------------------------

  if (!is.data.frame(data) && !is.matrix(data)) {
    stop(
      "`data` must be a data frame or matrix."
    )
  }

  data <- as.data.frame(
    data,
    check.names = FALSE
  )

  if (nrow(data) == 0) {
    stop(
      "`data` must contain at least one observation."
    )
  }

  if (ncol(data) == 0) {
    stop(
      "`data` must contain at least one trait."
    )
  }

  if (anyDuplicated(colnames(data))) {
    stop(
      "Trait names in `data` must be unique."
    )
  }


  # ---------------------------------------------------------------------------
  # 2. Calculate coverage for individual traits
  # ---------------------------------------------------------------------------

  n_observations <- nrow(
    data
  )

  trait_missing <- vapply(
    data,
    function(x) {
      sum(
        is.na(x)
      )
    },
    integer(1)
  )

  trait_available <-
    n_observations -
    trait_missing

  trait_coverage_pct <- round(
    trait_available /
      n_observations *
      100,
    2
  )


  trait_summary <- data.frame(
    Trait = names(
      trait_missing
    ),
    Available_count = trait_available,
    Missing_count = trait_missing,
    Trait_coverage_rate = paste0(
      trait_coverage_pct,
      " %"
    ),
    stringsAsFactors = FALSE,
    row.names = NULL
  )


  # ---------------------------------------------------------------------------
  # 3. Calculate complete-case coverage
  # ---------------------------------------------------------------------------

  complete_cases <- sum(
    stats::complete.cases(
      data
    )
  )

  complete_missing <-
    n_observations -
    complete_cases

  complete_coverage_pct <- round(
    complete_cases /
      n_observations *
      100,
    2
  )


  all_row <- data.frame(
    Trait = "All",
    Available_count = complete_cases,
    Missing_count = complete_missing,
    Trait_coverage_rate = paste0(
      complete_coverage_pct,
      " %"
    ),
    stringsAsFactors = FALSE,
    row.names = NULL
  )


  trait_summary <- rbind(
    trait_summary,
    all_row
  )


  # ---------------------------------------------------------------------------
  # 4. Calculate overall coverage
  # ---------------------------------------------------------------------------

  total_cells <-
    nrow(data) *
    ncol(data)

  total_missing <- sum(
    trait_missing
  )

  total_available <-
    total_cells -
    total_missing

  overall_coverage <- round(
    total_available /
      total_cells *
      100,
    2
  )


  message(
    "Overall trait coverage rate: ",
    overall_coverage,
    " %"
  )


  # ---------------------------------------------------------------------------
  # 5. Return results
  # ---------------------------------------------------------------------------

  trait_summary
}
