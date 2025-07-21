#' Calculate Trait Coverage Statistics
#'
#' This function calculates comprehensive coverage statistics for trait data,
#' including individual trait coverage rates, complete case coverage, and overall
#' data coverage. It provides both summary statistics and detailed breakdowns
#' of missing and available data.
#'
#' @param data A data frame containing trait data. Each column represents a trait
#'   and each row represents an observation (e.g., species, samples).
#'
#' @return A data frame with the following columns:
#'   \describe{
#'     \item{Trait}{Character. Names of traits plus an "All" row for complete cases}
#'     \item{Available_count}{Integer. Number of non-missing values for each trait}
#'     \item{Missing_count}{Integer. Number of missing (NA) values for each trait}
#'     \item{Trait_coverage_rate}{Character. Percentage of available data for each trait}
#'   }
#'   The "All" row shows statistics for complete cases (rows with no missing values).
#'
#' @details
#' The function performs the following calculations:
#' \itemize{
#'   \item \strong{Individual trait coverage}: For each trait, calculates the number
#'     and percentage of available (non-NA) values
#'   \item \strong{Complete case coverage}: Counts rows with no missing values across
#'     all traits and calculates the percentage
#'   \item \strong{Overall coverage}: Calculates the percentage of all cells in the
#'     dataset that contain non-missing values
#' }
#'
#' The function also prints the overall trait coverage rate to the console before
#' returning the detailed summary table.
#'
#' @examples
#' # Create sample trait data
#' trait_data <- data.frame(
#'   PlantHeight = c(1.2, 1.5, NA, 2.1, 1.8),
#'   LDMC = c(0.5, NA, 0.8, 1.2, 0.9),
#'   LA = c(15.2, 18.5, 12.3, NA, 16.7)
#' )
#'
#' # Calculate coverage statistics
#' coverage(trait_data)
#'
#' @export
coverage <- function(data) {
  trait_missing <- base::sapply(data, function(x) base::sum(base::is.na(x)))
  trait_coverage_pct <- base::round(((base::nrow(data) - trait_missing) / base::nrow(data)) * 100, 2)
  trait_summary <- base::data.frame(
    Trait = base::names(trait_missing),
    Available_count = base::nrow(data) - trait_missing,
    Missing_count = trait_missing,
    Trait_coverage_rate = base::paste0(trait_coverage_pct, " %"),
    row.names = NULL
  )
  complete_cases <- base::sum(stats::complete.cases(data))
  complete_missing <- base::nrow(data) - complete_cases
  complete_coverage_pct <- base::round((complete_cases / base::nrow(data)) * 100, 2)
  all_row <- base::data.frame(
    Trait = "All",
    Available_count = complete_cases,
    Missing_count = complete_missing,
    Trait_coverage_rate = base::paste0(complete_coverage_pct, " %"),
    row.names = NULL
  )
  trait_summary <- base::rbind(trait_summary, all_row)
  total_cells <- base::nrow(data) * base::ncol(data)
  total_missing <- base::sum(trait_missing)
  total_available <- total_cells - total_missing
  overall_coverage <- base::round((total_available / total_cells) * 100, 2)
  base::cat("Overall trait coverage rate:", base::paste0(overall_coverage, " %"), "\n\n")
  return(trait_summary)
}









