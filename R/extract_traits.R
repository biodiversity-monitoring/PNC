#' Extract Quantitative Traits for Specified Taxa
#'
#' `extract_traits()` retrieves selected quantitative traits for specified
#' species, genera, or families.
#'
#' Repeated records are first resolved separately for each species and trait.
#' For genus and family extraction, the resulting species-level values are
#' then aggregated so that each contributing species receives equal weight,
#' regardless of its number of records in the original dataset.
#'
#' @param sp.list A non-empty character vector containing unique taxon names
#'   to extract. Names must correspond to the taxonomic rank specified by
#'   `rank`.
#' @param dataset A non-empty data frame or tibble containing taxonomic
#'   information and quantitative trait data. A `"species"` column is required
#'   for all analyses. A `"genus"` or `"family"` column is additionally
#'   required when extracting traits at the corresponding rank.
#' @param rank Taxonomic rank used to match the requested taxa. Must be one of
#'   `"species"`, `"genus"`, or `"family"`. The default is `"species"`.
#' @param traits Character vector specifying the traits to extract. If `NULL`,
#'   all numeric columns other than `"species"`, `"genus"`, and `"family"`
#'   are extracted.
#' @param within_species Rule used when a species has multiple different
#'   finite records for the same trait. Must be `"error"`, `"mean"`, or
#'   `"median"`. The default, `"error"`, requires the user to resolve the
#'   records or explicitly select a summary rule.
#' @param among_species Rule used to aggregate species-level values within a
#'   genus or family. Must be `"mean"` or `"median"`. The default is `"mean"`.
#'   This argument is ignored when `rank = "species"`.
#'
#' @return A data frame with one row for each requested taxon, in the order
#'   supplied in `sp.list`, and one column for each retained numeric trait.
#'   Taxon names are stored as row names. Taxa not found in `dataset` are
#'   retained and returned with missing values for all traits.
#'
#' @details
#' Taxon names are matched exactly against the column corresponding to `rank`.
#' The function does not standardize taxonomic names, spelling, capitalization,
#' or synonyms.
#'
#' Only numeric traits are extracted. Requested traits that are absent or
#' non-numeric are ignored with a warning, and the function stops if no
#' requested numeric traits are available.
#'
#' For each requested trait, only records belonging to the requested taxa are
#' examined. Non-finite numeric values (`Inf`, `-Inf`, and `NaN`) in these
#' records are treated as missing and reported in a warning. Non-finite values
#' in unrequested taxa or traits do not affect the extraction.
#'
#' Within each species and trait, missing and non-finite records are excluded.
#' If no finite values remain, `NA` is returned. A single finite value, or
#' multiple identical finite values, is returned directly. Multiple different
#' finite values either trigger an error or are summarized using the mean or
#' median, according to `within_species`.
#'
#' For genus and family extraction, the resolved species-level values are
#' aggregated separately for each trait using `among_species`. Species without
#' a finite value for the focal trait are excluded from that trait's
#' aggregation. Each contributing species receives equal weight.
#'
#' The resulting genus or family value represents the central tendency of the
#' contributing species. When interpreting this aggregate, users should also
#' consider the number of contributing species and the variation among their
#' trait values.
#'
#' Each species must be assigned to only one requested genus or family.
#' Conflicting taxonomic assignments cause the function to stop. Values are
#' not inferred from higher taxonomic levels.
#'
#' Trait-specific validity constraints are not evaluated automatically. For
#' example, the function does not assume that zero is invalid because zero or
#' negative values may be valid for some quantitative traits.
#'
#' @examples
#' trait_data <- data.frame(
#'   species = c("Species a", "Species a", "Species b", "Species c"),
#'   genus = c("Genus1", "Genus1", "Genus1", "Genus2"),
#'   family = c("Family1", "Family1", "Family1", "Family2"),
#'   Trait1 = c(1, 3, 5, 2),
#'   Trait2 = c(4, 4, 6, NA)
#' )
#'
#' # Resolve repeated records using their mean
#' extract_traits(
#'   c("Species a", "Species c"),
#'   trait_data,
#'   traits = c("Trait1", "Trait2"),
#'   within_species = "mean"
#' )
#'
#' # Calculate genus-level means from resolved species-level values
#' extract_traits(
#'   c("Genus1", "Genus2"),
#'   trait_data,
#'   rank = "genus",
#'   traits = c("Trait1", "Trait2"),
#'   within_species = "mean",
#'   among_species = "mean"
#' )
#'
#' @export
extract_traits <- function(sp.list,
                           dataset,
                           rank = "species",
                           traits = NULL,
                           within_species = c(
                             "error",
                             "mean",
                             "median"
                           ),
                           among_species = c(
                             "mean",
                             "median"
                           )) {

  # Check inputs
  if (!is.character(sp.list) || length(sp.list) == 0L) {
    stop(
      "`sp.list` must be a non-empty character vector."
    )
  }

  if (anyNA(sp.list) || any(!nzchar(sp.list))) {
    stop(
      "`sp.list` must contain non-missing, non-empty taxon names."
    )
  }

  if (anyDuplicated(sp.list)) {
    stop(
      "Taxon names in `sp.list` must be unique."
    )
  }

  if (!is.data.frame(dataset) || nrow(dataset) == 0L) {
    stop(
      "`dataset` must be a non-empty data frame or tibble."
    )
  }

  if (anyDuplicated(names(dataset))) {
    stop(
      "Column names in `dataset` must be unique."
    )
  }

  rank <- match.arg(
    rank,
    choices = c(
      "species",
      "genus",
      "family"
    )
  )

  within_species <- match.arg(
    within_species
  )

  among_species <- match.arg(
    among_species
  )

  dataset <- as.data.frame(
    dataset,
    check.names = FALSE
  )

  required_taxonomy <- unique(
    c(
      "species",
      rank
    )
  )

  missing_taxonomy <- setdiff(
    required_taxonomy,
    names(dataset)
  )

  if (length(missing_taxonomy) > 0L) {
    stop(
      "`dataset` must contain the following taxonomic column(s): ",
      paste(
        missing_taxonomy,
        collapse = ", "
      ),
      "."
    )
  }

  taxonomy_cols <- intersect(
    c(
      "species",
      "genus",
      "family"
    ),
    names(dataset)
  )

  for (column in taxonomy_cols) {
    dataset[[column]] <- as.character(
      dataset[[column]]
    )
  }

  # Select numeric traits
  available_traits <- setdiff(
    names(dataset),
    taxonomy_cols
  )

  numeric_traits <- available_traits[
    vapply(
      dataset[available_traits],
      is.numeric,
      logical(1)
    )
  ]

  if (is.null(traits)) {

    trait_cols <- numeric_traits

    ignored_traits <- setdiff(
      available_traits,
      numeric_traits
    )

    if (length(ignored_traits) > 0L) {
      warning(
        "The following non-numeric columns were ignored: ",
        paste(
          ignored_traits,
          collapse = ", "
        ),
        ".",
        call. = FALSE
      )
    }

  } else {

    if (
      !is.character(traits) ||
      anyNA(traits) ||
      any(!nzchar(traits))
    ) {
      stop(
        "`traits` must be `NULL` or a character vector containing ",
        "non-missing, non-empty trait names."
      )
    }

    traits <- unique(
      traits
    )

    unavailable_traits <- setdiff(
      traits,
      numeric_traits
    )

    if (length(unavailable_traits) > 0L) {
      warning(
        "The following requested traits were absent or non-numeric and ",
        "were ignored: ",
        paste(
          unavailable_traits,
          collapse = ", "
        ),
        ".",
        call. = FALSE
      )
    }

    trait_cols <- traits[
      traits %in% numeric_traits
    ]
  }

  if (length(trait_cols) == 0L) {
    stop(
      "No numeric trait columns are available for extraction."
    )
  }

  # Match requested taxa
  matched_taxa <- sp.list[
    sp.list %in% dataset[[rank]]
  ]

  missing_taxa <- sp.list[
    !sp.list %in% dataset[[rank]]
  ]

  if (length(missing_taxa) > 0L) {
    message(
      length(missing_taxa),
      " ",
      rank,
      if (length(missing_taxa) == 1L) {
        " taxon"
      } else {
        " taxa"
      },
      " not found in the dataset: ",
      paste(
        missing_taxa,
        collapse = ", "
      )
    )
  }

  # Check taxonomic consistency for the requested genus or family
  if (rank != "species" && length(matched_taxa) > 0L) {

    requested_species <- unique(
      dataset$species[
        dataset[[rank]] %in% matched_taxa
      ]
    )

    valid_mapping <-
      dataset$species %in% requested_species &
      !is.na(dataset$species) &
      nzchar(dataset$species) &
      !is.na(dataset[[rank]]) &
      nzchar(dataset[[rank]])

    species_mapping <- unique(
      dataset[
        valid_mapping,
        c(
          "species",
          rank
        ),
        drop = FALSE
      ]
    )

    mapping_counts <- table(
      species_mapping$species
    )

    conflicting_species <- names(
      mapping_counts[
        mapping_counts > 1L
      ]
    )

    if (length(conflicting_species) > 0L) {
      stop(
        "The following species are assigned to more than one ",
        rank,
        ": ",
        paste(
          conflicting_species,
          collapse = ", "
        ),
        ". Resolve the taxonomic conflicts before extraction."
      )
    }
  }

  # Report non-finite values only in records used for this extraction
  matched_rows <- dataset[[rank]] %in% matched_taxa

  nonfinite_counts <- vapply(
    trait_cols,
    function(trait_name) {

      values <- dataset[[trait_name]][matched_rows]

      as.integer(
        sum(
          is.infinite(values) |
            is.nan(values)
        )
      )
    },
    integer(1)
  )

  nonfinite_counts <- nonfinite_counts[
    nonfinite_counts > 0L
  ]

  if (length(nonfinite_counts) > 0L) {

    nonfinite_summary <- paste0(
      names(nonfinite_counts),
      " (",
      nonfinite_counts,
      ")"
    )

    warning(
      "Non-finite values in matched records were treated as missing: ",
      paste(
        nonfinite_summary,
        collapse = ", "
      ),
      ".",
      call. = FALSE
    )
  }

  # Resolve repeated records within one species and trait
  resolve_species_value <- function(values,
                                    species_name,
                                    trait_name) {

    values <- values[
      is.finite(values)
    ]

    if (length(values) == 0L) {
      return(
        NA_real_
      )
    }

    distinct_values <- unique(
      values
    )

    if (length(distinct_values) == 1L) {
      return(
        as.numeric(distinct_values)
      )
    }

    if (within_species == "error") {
      stop(
        "Species `",
        species_name,
        "` has multiple different finite values for trait `",
        trait_name,
        "`. Set `within_species = \"mean\"` or `\"median\"`, ",
        "or resolve these records before extraction."
      )
    }

    if (within_species == "mean") {
      return(
        mean(values)
      )
    }

    stats::median(
      values
    )
  }

  # Initialize output
  result <- as.data.frame(
    matrix(
      NA_real_,
      nrow = length(sp.list),
      ncol = length(trait_cols),
      dimnames = list(
        sp.list,
        trait_cols
      )
    ),
    check.names = FALSE
  )

  # Extract species-level traits
  if (rank == "species") {

    for (species_name in matched_taxa) {

      species_rows <- which(
        dataset$species == species_name
      )

      for (trait_name in trait_cols) {

        result[
          species_name,
          trait_name
        ] <- resolve_species_value(
          dataset[[trait_name]][species_rows],
          species_name,
          trait_name
        )
      }
    }

  } else {

    # Aggregate species-level values within each genus or family
    for (taxon_name in matched_taxa) {

      taxon_rows <- which(
        dataset[[rank]] == taxon_name
      )

      species_names <- unique(
        dataset$species[taxon_rows]
      )

      species_names <- species_names[
        !is.na(species_names) &
          nzchar(species_names)
      ]

      for (trait_name in trait_cols) {

        species_values <- rep(
          NA_real_,
          length(species_names)
        )

        for (i in seq_along(species_names)) {

          species_name <- species_names[i]

          species_rows <- taxon_rows[
            which(
              dataset$species[taxon_rows] ==
                species_name
            )
          ]

          species_values[i] <- resolve_species_value(
            dataset[[trait_name]][species_rows],
            species_name,
            trait_name
          )
        }

        species_values <- species_values[
          is.finite(species_values)
        ]

        if (length(species_values) == 0L) {
          next
        }

        if (among_species == "mean") {

          result[
            taxon_name,
            trait_name
          ] <- mean(
            species_values
          )

        } else {

          result[
            taxon_name,
            trait_name
          ] <- stats::median(
            species_values
          )
        }
      }
    }
  }

  # Ensure that numeric output contains only finite values or NA
  result[] <- lapply(
    result,
    function(values) {

      values[
        !is.finite(values)
      ] <- NA_real_

      values
    }
  )

  result
}
