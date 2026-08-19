#' Calculate Trait Coverage Statistics
#'
#' `coverage()` summarizes the completeness of a trait dataset by calculating
#' coverage for individual traits, complete cases across all traits, and the
#' overall proportion of non-missing values.
#'
#' @param data A data frame or matrix containing trait data. Rows represent
#'   observations such as species, and columns represent traits.
#'
#' @return A data frame with the following columns:
#'   \describe{
#'     \item{Trait}{Trait name, with an additional `"All"` row representing
#'       complete cases across all traits.}
#'     \item{Available_count}{Number of non-missing values.}
#'     \item{Missing_count}{Number of missing values.}
#'     \item{Trait_coverage_rate}{Percentage of available values, formatted
#'       as a character string.}
#'   }
#'
#' @details
#' Coverage is calculated in three ways:
#' \itemize{
#'   \item Individual trait coverage: percentage of observations with a
#'     non-missing value for each trait.
#'   \item Complete-case coverage: percentage of observations with non-missing
#'     values for all traits, reported in the `"All"` row.
#'   \item Overall coverage: percentage of all cells in the dataset that are
#'     non-missing. This value is reported as a message when the function runs.
#' }
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
