#' Merge Two Trait Datasets by Species
#'
#' `merge_dataset()` combines two trait datasets using the `species` column as
#' the matching key. Species and variables occurring in either dataset are
#' retained, and missing or conflicting values are resolved according to the
#' selected priority rule.
#'
#' @param main_data A data frame containing the primary dataset. It must include
#'   a column named `"species"`.
#' @param additional_data A data frame containing the additional dataset. It
#'   must include a column named `"species"`.
#' @param priority Character string specifying how conflicts are resolved when
#'   both datasets contain non-missing values for the same species and variable.
#'   Available options are:
#'   \itemize{
#'     \item `"main"`: retain the value from `main_data` (default).
#'     \item `"additional"`: retain the value from `additional_data`.
#'     \item `"mean"`: use the arithmetic mean when the corresponding columns
#'       in both datasets are numeric; otherwise retain the value from
#'       `main_data`.
#'   }
#'
#' @return A data frame containing all unique species and all variables present
#'   in either input dataset. The `species` column is placed first. Other
#'   columns retain their order from `main_data`, followed by columns occurring
#'   only in `additional_data`.
#'
#' @details
#' For a given species and variable, a non-missing value is always retained when
#' the corresponding value in the other dataset is missing. The `priority`
#' argument is used only when both datasets contain non-missing values.
#'
#' If a species occurs more than once within an input dataset, only its first
#' occurrence is used.
#'
#' With `priority = "mean"`, averaging is performed only when the corresponding
#' columns in both datasets are numeric.
#'
#' @examples
#' main_data <- data.frame(
#'   species = c(
#'     "Abies alba",
#'     "Coussapoa trinervia",
#'     "Crataegus monogyna"
#'   ),
#'   genus = c(
#'     "Abies",
#'     "Coussapoa",
#'     "Crataegus"
#'   ),
#'   family = c(
#'     "Pinaceae",
#'     "Urticaceae",
#'     "Rosaceae"
#'   ),
#'   LA = c(NA, 2050.24, 449.15),
#'   LeafN = c(13.10, 14.52, 17.46),
#'   SeedMass = c(53.64, NA, 95.92),
#'   stringsAsFactors = FALSE
#' )
#'
#' additional_data <- data.frame(
#'   species = c(
#'     "Abies alba",
#'     "Corydalis solida"
#'   ),
#'   genus = c(
#'     "Abies",
#'     "Corydalis"
#'   ),
#'   family = c(
#'     "Pinaceae",
#'     "Papaveraceae"
#'   ),
#'   LA = c(25.58, NA),
#'   LMA = c(0.19, 0.20),
#'   PlantHeight = c(53.66, 0.14),
#'   stringsAsFactors = FALSE
#' )
#'
#' merge_dataset(
#'   main_data,
#'   additional_data
#' )
#'
#' merge_dataset(
#'   main_data,
#'   additional_data,
#'   priority = "mean"
#' )
#'
#' @export
merge_dataset <- function(main_data,
                          additional_data,
                          priority = "main") {

  # ---------------------------------------------------------------------------
  # 1. Check inputs
  # ---------------------------------------------------------------------------

  if (!is.data.frame(main_data)) {
    stop("`main_data` must be a data frame.")
  }

  if (!is.data.frame(additional_data)) {
    stop("`additional_data` must be a data frame.")
  }

  if (!"species" %in% names(main_data)) {
    stop("`main_data` must contain a `species` column.")
  }

  if (!"species" %in% names(additional_data)) {
    stop("`additional_data` must contain a `species` column.")
  }

  priority <- match.arg(
    priority,
    choices = c("main", "additional", "mean")
  )

  if (
    anyNA(main_data$species) ||
    any(!nzchar(as.character(main_data$species)))
  ) {
    stop(
      "`main_data$species` must contain non-missing, non-empty species names."
    )
  }

  if (
    anyNA(additional_data$species) ||
    any(!nzchar(as.character(additional_data$species)))
  ) {
    stop(
      "`additional_data$species` must contain non-missing, non-empty species names."
    )
  }


  # ---------------------------------------------------------------------------
  # 2. Prepare datasets
  # ---------------------------------------------------------------------------

  main_data$species <- as.character(main_data$species)
  additional_data$species <- as.character(additional_data$species)

  if (anyDuplicated(main_data$species)) {

    warning(
      "Duplicate species were found in `main_data`; ",
      "only the first occurrence of each species was used.",
      call. = FALSE
    )

    main_data <- main_data[
      !duplicated(main_data$species),
      ,
      drop = FALSE
    ]
  }

  if (anyDuplicated(additional_data$species)) {

    warning(
      "Duplicate species were found in `additional_data`; ",
      "only the first occurrence of each species was used.",
      call. = FALSE
    )

    additional_data <- additional_data[
      !duplicated(additional_data$species),
      ,
      drop = FALSE
    ]
  }

  main_data[] <- lapply(
    main_data,
    function(x) {
      if (is.factor(x)) {
        as.character(x)
      } else {
        x
      }
    }
  )

  additional_data[] <- lapply(
    additional_data,
    function(x) {
      if (is.factor(x)) {
        as.character(x)
      } else {
        x
      }
    }
  )


  # ---------------------------------------------------------------------------
  # 3. Define species and variables
  # ---------------------------------------------------------------------------

  all_species <- unique(
    c(
      main_data$species,
      additional_data$species
    )
  )

  all_columns <- unique(
    c(
      names(main_data),
      names(additional_data)
    )
  )

  trait_columns <- setdiff(
    all_columns,
    "species"
  )

  main_index <- match(
    all_species,
    main_data$species
  )

  additional_index <- match(
    all_species,
    additional_data$species
  )

  result <- data.frame(
    species = all_species,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )


  # ---------------------------------------------------------------------------
  # 4. Merge variables
  # ---------------------------------------------------------------------------

  for (column in trait_columns) {

    in_main <- column %in% names(main_data)
    in_additional <- column %in% names(additional_data)


    # Variable occurs only in main_data
    if (in_main && !in_additional) {

      result[[column]] <- main_data[[column]][main_index]

      next
    }


    # Variable occurs only in additional_data
    if (!in_main && in_additional) {

      result[[column]] <- additional_data[[column]][additional_index]

      next
    }


    # Variable occurs in both datasets
    main_values <- main_data[[column]][main_index]

    additional_values <- additional_data[[column]][additional_index]

    main_available <- !is.na(main_values)

    additional_available <- !is.na(additional_values)


    # -------------------------------------------------------------------------
    # Mean priority for numeric columns
    # -------------------------------------------------------------------------

    if (
      priority == "mean" &&
      is.numeric(main_data[[column]]) &&
      is.numeric(additional_data[[column]])
    ) {

      merged_values <- rep(
        NA_real_,
        length(all_species)
      )

      use_main <- main_available &
        !additional_available

      merged_values[use_main] <- main_values[use_main]

      use_additional <- !main_available &
        additional_available

      merged_values[use_additional] <-
        additional_values[use_additional]

      use_mean <- main_available &
        additional_available

      merged_values[use_mean] <-
        (
          main_values[use_mean] +
            additional_values[use_mean]
        ) / 2

      result[[column]] <- merged_values

      next
    }


    # -------------------------------------------------------------------------
    # Main priority
    # -------------------------------------------------------------------------

    if (
      priority == "main" ||
      priority == "mean"
    ) {

      merged_values <- main_values

      use_additional <- !main_available &
        additional_available

      merged_values[use_additional] <-
        additional_values[use_additional]

      result[[column]] <- merged_values

      next
    }


    # -------------------------------------------------------------------------
    # Additional priority
    # -------------------------------------------------------------------------

    merged_values <- additional_values

    use_main <- !additional_available &
      main_available

    merged_values[use_main] <-
      main_values[use_main]

    result[[column]] <- merged_values
  }


  # ---------------------------------------------------------------------------
  # 5. Return results
  # ---------------------------------------------------------------------------

  result
}
