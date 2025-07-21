#' Merge Two Datasets Based on Species Column
#'
#' This function merges two data frames based on the 'species' column, handling
#' missing values and column differences intelligently. It provides flexible
#' options for resolving conflicts when the same species appears in both datasets.
#'
#' @param main_data A data frame containing the primary dataset. Must include a 'species' column.
#' @param additional_data A data frame containing the secondary dataset. Must include a 'species' column.
#' @param priority A character string specifying how to handle conflicts when both datasets
#'   contain non-missing values for the same species and column. Options are:
#'   \itemize{
#'     \item \code{"main"} (default): Use values from main_data
#'     \item \code{"additional"}: Use values from additional_data
#'     \item \code{"mean"}: Calculate mean for numeric values, use main_data for non-numeric
#'   }
#'
#' @return A data frame containing all unique species from both input datasets, with
#'   all columns from both datasets. The 'species' column is placed first, followed
#'   by all other columns in alphabetical order.
#'
#' @details
#' The function performs the following operations:
#' \itemize{
#'   \item Combines all unique species from both datasets
#'   \item Includes all columns from both datasets
#'   \item Handles missing values by using available non-missing values
#'   \item Resolves conflicts based on the specified priority
#'   \item For duplicate species within a dataset, only the first occurrence is used
#' }
#'
#' @examples
#' # Create sample datasets
#' main_data <- data.frame(
#'   species = c("Abies alba", "Coussapoa trinervia", "Crataegus monogyna"),
#'   genus = c("Abies", "Coussapoa", "Crataegus"),
#'   family = c("Pinaceae", "Urticaceae", "Rosaceae"),
#'   LA = c(NA, 2050.24, 449.15),
#'   LeafN = c(13.10, 14.52, 17.46),
#'   Seedmass = c(53.64, NA, 95.92),
#'   stringsAsFactors = FALSE
#' )
#'
#' additional_data <- data.frame(
#'   species = c("Abies alba", "Corydalis solida"),
#'   genus = c("Abies", "Corydalis"),
#'   family = c("Pinaceae", "Papaveraceae"),
#'   LA = c(25.58, NA),
#'   LMA = c(0.19, 0.2),
#'   PlantHeight = c(53.66, 0.14),
#'   stringsAsFactors = FALSE
#' )
#'
#' # Merge with main data priority (default)
#' merge_dataset(main_data, additional_data)
#'
#' @export
#'
#' @importFrom stats na.omit
#'
#' @note
#' \itemize{
#'   \item Both input datasets must contain a 'species' column
#'   \item If a species appears multiple times in a dataset, only the first occurrence is used
#'   \item When priority is "mean", non-numeric values default to main_data values
#'   \item The function preserves the original data types of columns
#' }
merge_dataset <- function(main_data, additional_data, priority = "main") {
  if (!is.data.frame(main_data) || !is.data.frame(additional_data)) {
    stop("main_data and additional_data must be data frames.", call. = FALSE)
  }
  if (!"species" %in% names(main_data) || !"species" %in% names(additional_data)) {
    stop("Both datasets must contain the 'species' column.", call. = FALSE)
  }
  if (!priority %in% c("main", "additional", "mean")) {
    stop("The priority parameter must be one of 'main', 'additional', or 'mean'.",
         call. = FALSE)
  }
  all_species <- base::unique(c(main_data$species, additional_data$species))
  all_columns <- base::unique(c(names(main_data), names(additional_data)))
  result <- data.frame(
    matrix(NA, nrow = length(all_species), ncol = length(all_columns)),
    stringsAsFactors = FALSE
  )
  names(result) <- all_columns
  result$species <- all_species
  merge_values <- function(main_val, add_val, priority) {
    if (base::is.na(main_val) && base::is.na(add_val)) {
      return(NA)
    }
    if (base::is.na(main_val) && !base::is.na(add_val)) {
      return(add_val)
    }
    if (!base::is.na(main_val) && base::is.na(add_val)) {
      return(main_val)
    }
    if (!base::is.na(main_val) && !base::is.na(add_val)) {
      switch(priority,
             "main" = main_val,
             "additional" = add_val,
             "mean" = {
               if (is.numeric(main_val) && is.numeric(add_val)) {
                 (main_val + add_val) / 2
               } else {
                 main_val
               }
             }
      )
    }
  }
  for (species in all_species) {
    main_rows <- base::which(main_data$species == species)
    add_rows <- base::which(additional_data$species == species)
    main_row <- if (length(main_rows) > 0) main_rows[1] else NA
    add_row <- if (length(add_rows) > 0) add_rows[1] else NA
    result_row <- base::which(result$species == species)
    for (col in all_columns) {
      if (col == "species") {
        next
      }
      main_val <- if (!base::is.na(main_row) && col %in% names(main_data)) {
        main_data[main_row, col]
      } else {
        NA
      }
      add_val <- if (!base::is.na(add_row) && col %in% names(additional_data)) {
        additional_data[add_row, col]
      } else {
        NA
      }
      result[result_row, col] <- merge_values(main_val, add_val, priority)
    }
  }
  species_col <- result$species
  other_cols <- result[, !names(result) %in% "species", drop = FALSE]
  result <- data.frame(
    species = species_col,
    other_cols,
    stringsAsFactors = FALSE
  )
  return(result)
}

