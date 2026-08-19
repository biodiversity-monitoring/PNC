#' Extract Traits from a Trait Database
#'
#' `extract_traits()` extracts trait data for a specified list of taxa from a
#' trait database at the species, genus, or family level.
#'
#' At the species level, trait values are taken from the first matching record.
#' At the genus and family levels, numeric traits are averaged across all
#' available records within each taxon. For non-numeric traits, the first
#' non-missing value is returned.
#'
#' @param sp.list Character vector containing the taxa to extract. Taxon names
#'   should correspond to the taxonomic rank specified by `rank`.
#' @param dataset A data frame or tibble containing taxonomic information and
#'   trait data. The column corresponding to `rank` must be present. Columns
#'   named `"species"`, `"genus"`, and `"family"` are treated as taxonomic
#'   columns rather than traits.
#' @param rank Taxonomic rank used for matching taxa. Must be one of
#'   `"species"`, `"genus"`, or `"family"`. Default is `"species"`.
#' @param traits Character vector specifying the traits to extract. If `NULL`,
#'   all non-taxonomic columns are extracted.
#'
#' @return A data frame with taxa as row names and traits as columns. Taxa not
#'   found in the dataset are retained as rows with missing trait values.
#'
#' @details
#' At the species level, the first matching record is used when duplicate
#' species records occur in the dataset.
#'
#' At the genus and family levels, numeric traits are averaged across available
#' records using `na.rm = TRUE`. If all values for a numeric trait are missing,
#' `NA` is returned. For non-numeric traits, the first non-missing value is
#' returned.
#'
#' Requested traits that are not present in the dataset are ignored with a
#' warning. The function stops if none of the requested traits are available.
#'
#' @examples
#' data(TRY)
#'
#' species_list <- c(
#'   "Acaena novae-zelandiae",
#'   "Adiantum capillus-veneris",
#'   "Zuelania guidonia"
#' )
#'
#' # Extract all available traits
#' extract_traits(
#'   species_list,
#'   TRY,
#'   rank = "species"
#' )
#'
#' # Extract selected traits
#' extract_traits(
#'   species_list,
#'   TRY,
#'   rank = "species",
#'   traits = c(
#'     "LA",
#'     "LMA",
#'     "LeafN",
#'     "PlantHeight",
#'     "SeedMass",
#'     "SSD"
#'   )
#' )
#'
#' # Extract traits at the genus level
#' genus_list <- c(
#'   "Acaena",
#'   "Adiantum"
#' )
#'
#' extract_traits(
#'   genus_list,
#'   TRY,
#'   rank = "genus",
#'   traits = c(
#'     "LDMC",
#'     "PlantHeight",
#'     "SeedMass"
#'   )
#' )
#'
#' @export
extract_traits <- function(sp.list,
                           dataset,
                           rank = "species",
                           traits = NULL) {

  # ---------------------------------------------------------------------------
  # 1. Check inputs
  # ---------------------------------------------------------------------------

  if (!is.character(sp.list)) {
    stop(
      "`sp.list` must be a character vector."
    )
  }

  if (length(sp.list) == 0) {
    stop(
      "`sp.list` must contain at least one taxon."
    )
  }

  if (
    anyNA(sp.list) ||
    any(!nzchar(sp.list))
  ) {
    stop(
      "`sp.list` must contain non-missing, non-empty taxon names."
    )
  }

  if (anyDuplicated(sp.list)) {
    stop(
      "Taxon names in `sp.list` must be unique."
    )
  }

  if (!is.data.frame(dataset)) {
    stop(
      "`dataset` must be a data frame or tibble."
    )
  }

  if (nrow(dataset) == 0 || ncol(dataset) == 0) {
    stop(
      "`dataset` must contain data."
    )
  }

  # Convert tibbles and other data-frame subclasses to a standard data.frame.
  # This allows stable row-name handling in the returned object.
  dataset <- as.data.frame(
    dataset,
    check.names = FALSE
  )

  rank <- match.arg(
    rank,
    choices = c(
      "species",
      "genus",
      "family"
    )
  )

  if (!rank %in% names(dataset)) {
    stop(
      "`dataset` must contain a `",
      rank,
      "` column."
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

  available_traits <- setdiff(
    names(dataset),
    taxonomy_cols
  )

  if (length(available_traits) == 0) {
    stop(
      "`dataset` contains no trait columns."
    )
  }


  # ---------------------------------------------------------------------------
  # 2. Select traits
  # ---------------------------------------------------------------------------

  if (is.null(traits)) {

    trait_cols <- available_traits

  } else {

    if (!is.character(traits)) {
      stop(
        "`traits` must be a character vector or `NULL`."
      )
    }

    if (
      anyNA(traits) ||
      any(!nzchar(traits))
    ) {
      stop(
        "`traits` must contain non-missing, non-empty trait names."
      )
    }

    traits <- unique(
      traits
    )

    invalid_traits <- setdiff(
      traits,
      available_traits
    )

    if (length(invalid_traits) > 0) {
      warning(
        "The following requested traits were not found and were ignored: ",
        paste(
          invalid_traits,
          collapse = ", "
        ),
        ".",
        call. = FALSE
      )
    }

    trait_cols <- intersect(
      traits,
      available_traits
    )

    if (length(trait_cols) == 0) {
      stop(
        "None of the requested traits are available in `dataset`."
      )
    }
  }


  # ---------------------------------------------------------------------------
  # 3. Initialize output
  # ---------------------------------------------------------------------------

  # Start from one existing row so that the original column classes are
  # retained, then replace all values with missing values.
  result <- dataset[
    rep(
      1L,
      length(sp.list)
    ),
    trait_cols,
    drop = FALSE
  ]

  for (trait_name in trait_cols) {

    if (is.numeric(result[[trait_name]])) {

      result[[trait_name]] <- rep(
        NA_real_,
        length(sp.list)
      )

    } else if (is.integer(result[[trait_name]])) {

      result[[trait_name]] <- rep(
        NA_integer_,
        length(sp.list)
      )

    } else if (is.logical(result[[trait_name]])) {

      result[[trait_name]] <- rep(
        NA,
        length(sp.list)
      )

    } else {

      result[[trait_name]] <- rep(
        NA_character_,
        length(sp.list)
      )
    }
  }

  rownames(result) <- sp.list


  # ---------------------------------------------------------------------------
  # 4. Match taxa
  # ---------------------------------------------------------------------------

  trait_taxa <- dataset[[rank]]

  matched_taxa <- intersect(
    sp.list,
    trait_taxa
  )

  missing_taxa <- setdiff(
    sp.list,
    trait_taxa
  )

  if (length(missing_taxa) > 0) {
    message(
      length(missing_taxa),
      " ",
      rank,
      if (length(missing_taxa) == 1) {
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


  # ---------------------------------------------------------------------------
  # 5. Extract traits
  # ---------------------------------------------------------------------------

  for (taxon in matched_taxa) {

    taxon_rows <- which(
      dataset[[rank]] == taxon
    )


    # -------------------------------------------------------------------------
    # Species level
    # -------------------------------------------------------------------------

    if (rank == "species") {

      result[
        taxon,
        trait_cols
      ] <- dataset[
        taxon_rows[1],
        trait_cols,
        drop = FALSE
      ]

      next
    }


    # -------------------------------------------------------------------------
    # Genus or family level
    # -------------------------------------------------------------------------

    taxon_data <- dataset[
      taxon_rows,
      trait_cols,
      drop = FALSE
    ]

    for (trait_name in trait_cols) {

      values <- taxon_data[[trait_name]]


      # Numeric traits: arithmetic mean across available records
      if (is.numeric(values)) {

        if (all(is.na(values))) {

          result[
            taxon,
            trait_name
          ] <- NA_real_

        } else {

          result[
            taxon,
            trait_name
          ] <- mean(
            values,
            na.rm = TRUE
          )
        }

        next
      }


      # Non-numeric traits: first available value
      non_missing <- which(
        !is.na(values)
      )

      if (length(non_missing) > 0) {

        result[
          taxon,
          trait_name
        ] <- values[
          non_missing[1]
        ]
      }
    }
  }


  # ---------------------------------------------------------------------------
  # 6. Replace NaN with NA
  # ---------------------------------------------------------------------------

  result[] <- lapply(
    result,
    function(x) {

      if (is.numeric(x)) {
        x[is.nan(x)] <- NA_real_
      }

      x
    }
  )


  # ---------------------------------------------------------------------------
  # 7. Return results
  # ---------------------------------------------------------------------------

  result
}
