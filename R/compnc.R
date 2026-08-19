#' Estimate Phylogenetic Signal Across Multiple Communities
#'
#' `compnc()` estimates phylogenetic signal in functional traits across
#' multiple ecological communities. The function identifies species present
#' in each community, matches them with trait data and a phylogenetic tree,
#' and estimates phylogenetic signal using Pagel's lambda or Blomberg's K.
#'
#' When PCA axes are requested, a single PCA is performed using the pooled
#' species set represented across the communities, trait data, and phylogeny.
#' The resulting common PCA space is then used for all communities, ensuring
#' that a given PCA axis represents the same multivariate trait dimension
#' across communities.
#'
#' Phylogenetic signal describes the tendency for related species to resemble
#' one another in trait values. It can provide evidence relevant to
#' phylogenetic niche conservatism, but does not by itself constitute a
#' direct test of phylogenetic niche conservatism.
#'
#' @param com A community matrix or data frame with communities as rows and
#'   species as columns. Cell values represent species abundance or occurrence.
#' @param trait_data A data frame or matrix containing continuous trait data,
#'   with species as rows and traits as columns. Species names must be stored
#'   as row names.
#' @param phylo_tree A phylogenetic tree object of class `"phylo"`.
#' @param methods Character vector specifying the phylogenetic signal metrics
#'   to calculate. Available options are `"lambda"` and `"K"`.
#' @param pca_axes Character vector specifying PCA axes to include, for example
#'   `c("PC1", "PC2")`. Set to `NULL` to skip PCA.
#' @param sig_levels Numeric vector of three increasing significance thresholds
#'   used to assign significance symbols. Default is
#'   `c(0.001, 0.01, 0.05)`.
#' @param min_abundance Minimum abundance required for a species to be
#'   considered present in a community. A species is considered present when
#'   its abundance is greater than `min_abundance`. Default is 0.
#' @param nsim Number of randomizations used for significance testing of
#'   Blomberg's K. Default is 1000.
#' @param verbose Logical. If `TRUE`, warnings from PCA or phylogenetic signal
#'   estimation are reported.
#'
#' @return A data frame with one row for each
#'   community-trait-method combination:
#'   \itemize{
#'     \item `plot`: community name.
#'     \item `trait`: trait name or PCA axis.
#'     \item `coverage`: percentage of species present in the community with
#'       available values for the analyzed trait or PCA axis.
#'     \item `n_sp`: number of species actually included in the phylogenetic
#'       signal analysis after matching trait data with the phylogeny.
#'     \item `signal`: estimated phylogenetic signal.
#'     \item `p`: P value.
#'     \item `significance`: significance symbol.
#'     \item `method`: phylogenetic signal metric.
#'     \item `n_sp_in_plot`: total number of species present in the community
#'       according to `min_abundance`.
#'   }
#'
#' @details
#' Original traits are analyzed independently, so missing data in one trait
#' do not affect the species included for another trait.
#'
#' When PCA axes are requested, the PCA is based on species that occur in at
#' least one community, are represented in the phylogeny, are present in
#' `trait_data`, and have complete data across all supplied traits. The same
#' PCA model is then used for every community.
#'
#' Phylogenetic signal is not estimated for a community-trait combination when
#' fewer than four species have both trait values and phylogenetic information.
#' This threshold is used as a practical computational safeguard rather than
#' as a universal minimum sample-size requirement for phylogenetic signal
#' estimation.
#'
#' @examples
#' \donttest{
#' data(HimalayanBirds)
#' data(AVONET)
#'
#' sp <- colnames(HimalayanBirds$com)
#'
#' subtraits <- extract_traits(
#'   sp,
#'   AVONET,
#'   rank = "species"
#' )
#'
#' # Pagel's lambda with a common PCA space
#' compnc(
#'   HimalayanBirds$com,
#'   subtraits,
#'   HimalayanBirds$phy_species,
#'   methods = "lambda"
#' )
#'
#' # Pagel's lambda without PCA
#' compnc(
#'   HimalayanBirds$com,
#'   subtraits,
#'   HimalayanBirds$phy_species,
#'   methods = "lambda",
#'   pca_axes = NULL
#' )
#' }
#'
#' @references
#' Münkemüller, T., Lavergne, S., Bzeznik, B., Dray, S., Jombart, T.,
#' Schiffers, K. and Thuiller, W. (2012). How to measure and test
#' phylogenetic signal. Methods in Ecology and Evolution, 3, 743-756.
#' \doi{10.1111/j.2041-210X.2012.00196.x}
#'
#' @export
compnc <- function(com,
                   trait_data,
                   phylo_tree,
                   methods = "lambda",
                   pca_axes = c("PC1", "PC2"),
                   sig_levels = c(0.001, 0.01, 0.05),
                   min_abundance = 0,
                   nsim = 1000,
                   verbose = TRUE) {

  # ---------------------------------------------------------------------------
  # 1. Check inputs
  # ---------------------------------------------------------------------------

  if (!is.data.frame(com) && !is.matrix(com)) {
    stop("`com` must be a data frame or matrix.")
  }

  if (!is.data.frame(trait_data) && !is.matrix(trait_data)) {
    stop("`trait_data` must be a data frame or matrix.")
  }

  if (!inherits(phylo_tree, "phylo")) {
    stop("`phylo_tree` must be an object of class `phylo`.")
  }


  # Preserve names before conversion to data frames
  species_names <- colnames(com)
  plot_names <- rownames(com)
  trait_species <- rownames(trait_data)


  if (is.null(species_names)) {
    stop(
      "Species names must be provided as column names of `com`."
    )
  }

  if (
    is.null(trait_species) ||
    identical(
      trait_species,
      as.character(seq_len(nrow(trait_data)))
    )
  ) {
    stop(
      "`trait_data` must have species names as row names."
    )
  }


  com <- as.data.frame(
    com,
    check.names = FALSE
  )

  trait_data <- as.data.frame(
    trait_data,
    check.names = FALSE
  )


  if (nrow(com) == 0 || ncol(com) == 0) {
    stop(
      "`com` must contain at least one community and one species."
    )
  }

  if (nrow(trait_data) == 0 || ncol(trait_data) == 0) {
    stop(
      "`trait_data` must contain at least one species and one trait."
    )
  }


  # Generate community names when none were supplied
  if (
    is.null(plot_names) ||
    identical(
      plot_names,
      as.character(seq_len(nrow(com)))
    )
  ) {

    plot_names <- paste0(
      "plot",
      sprintf(
        "%03d",
        seq_len(nrow(com))
      )
    )

    rownames(com) <- plot_names

  } else {

    if (anyDuplicated(plot_names)) {
      stop(
        "Community names in `com` must be unique."
      )
    }

    rownames(com) <- plot_names
  }


  if (anyDuplicated(species_names)) {
    stop(
      "Species names in `com` must be unique."
    )
  }

  if (anyDuplicated(trait_species)) {
    stop(
      "Species names in `trait_data` must be unique."
    )
  }

  if (anyDuplicated(colnames(trait_data))) {
    stop(
      "Trait names in `trait_data` must be unique."
    )
  }

  if (anyDuplicated(phylo_tree$tip.label)) {
    stop(
      "Tip labels in `phylo_tree` must be unique."
    )
  }


  numeric_com <- vapply(
    com,
    is.numeric,
    logical(1)
  )

  if (!all(numeric_com)) {
    stop(
      "All columns in `com` must be numeric."
    )
  }


  numeric_traits <- vapply(
    trait_data,
    is.numeric,
    logical(1)
  )

  if (!all(numeric_traits)) {
    stop(
      "All trait columns must be numeric. Non-numeric traits: ",
      paste(
        names(numeric_traits)[!numeric_traits],
        collapse = ", "
      )
    )
  }


  # Check non-finite community values
  non_finite_com <- vapply(
    com,
    function(x) {
      any(
        !is.na(x) &
          !is.finite(x)
      )
    },
    logical(1)
  )

  if (any(non_finite_com)) {
    stop(
      "`com` contains non-finite abundance values."
    )
  }


  # Community abundance or occurrence values should not be negative
  negative_com <- vapply(
    com,
    function(x) {
      any(
        !is.na(x) &
          x < 0
      )
    },
    logical(1)
  )

  if (any(negative_com)) {
    stop(
      "`com` contains negative abundance or occurrence values."
    )
  }


  # Check non-finite trait values
  non_finite_traits <- vapply(
    trait_data,
    function(x) {
      any(
        !is.na(x) &
          !is.finite(x)
      )
    },
    logical(1)
  )

  if (any(non_finite_traits)) {
    stop(
      "Trait data contain non-finite values in: ",
      paste(
        names(non_finite_traits)[non_finite_traits],
        collapse = ", "
      ),
      "."
    )
  }


  if (
    length(min_abundance) != 1 ||
    !is.finite(min_abundance) ||
    min_abundance < 0
  ) {
    stop(
      "`min_abundance` must be a single non-negative finite number."
    )
  }


  methods <- match.arg(
    methods,
    choices = c("lambda", "K"),
    several.ok = TRUE
  )


  if (
    !is.null(pca_axes) &&
    !is.character(pca_axes)
  ) {
    stop(
      "`pca_axes` must be a character vector or `NULL`."
    )
  }


  if (!is.null(pca_axes)) {

    if (
      anyNA(pca_axes) ||
      any(!nzchar(pca_axes))
    ) {
      stop(
        "`pca_axes` must contain non-missing, non-empty names."
      )
    }

    pca_axes <- unique(
      pca_axes
    )

    duplicated_names <- intersect(
      pca_axes,
      colnames(trait_data)
    )

    if (length(duplicated_names) > 0) {
      stop(
        "Requested PCA axis names already exist in `trait_data`: ",
        paste(
          duplicated_names,
          collapse = ", "
        ),
        "."
      )
    }
  }


  if (
    length(sig_levels) != 3 ||
    any(!is.finite(sig_levels)) ||
    any(sig_levels <= 0) ||
    any(sig_levels >= 1) ||
    any(diff(sig_levels) <= 0)
  ) {
    stop(
      "`sig_levels` must contain three increasing values between 0 and 1."
    )
  }


  if (
    length(nsim) != 1 ||
    !is.finite(nsim) ||
    nsim < 1 ||
    nsim %% 1 != 0
  ) {
    stop(
      "`nsim` must be a positive integer."
    )
  }


  if (
    length(verbose) != 1 ||
    !is.logical(verbose) ||
    is.na(verbose)
  ) {
    stop(
      "`verbose` must be TRUE or FALSE."
    )
  }


  nsim <- as.integer(
    nsim
  )


  # ---------------------------------------------------------------------------
  # 2. Define study species and common PCA space
  # ---------------------------------------------------------------------------

  # Species occurring in at least one community
  species_present <- vapply(
    com,
    function(x) {
      any(
        !is.na(x) &
          x > min_abundance
      )
    },
    logical(1)
  )

  study_species <- colnames(com)[
    species_present
  ]


  if (length(study_species) == 0) {
    stop(
      "No species are present in `com` at the specified `min_abundance`."
    )
  }


  # Species eligible to define the common PCA space
  reference_species <- Reduce(
    intersect,
    list(
      study_species,
      rownames(trait_data),
      phylo_tree$tip.label
    )
  )


  if (length(reference_species) == 0) {
    stop(
      "No species are shared among `com`, `trait_data`, and `phylo_tree`."
    )
  }


  analysis_data <- trait_data

  pca_model <- NULL
  pca_results <- NULL
  pca_failed <- FALSE


  if (length(pca_axes) > 0) {

    pca_data <- trait_data[
      reference_species,
      ,
      drop = FALSE
    ]

    complete_rows <- stats::complete.cases(
      pca_data
    )

    complete_data <- pca_data[
      complete_rows,
      ,
      drop = FALSE
    ]


    if (nrow(complete_data) < 4) {

      pca_failed <- TRUE

      if (verbose) {
        warning(
          "PCA was skipped because fewer than four species ",
          "shared among `com`, `trait_data`, and `phylo_tree` ",
          "had complete trait data."
        )
      }

    } else {

      pca_model <- tryCatch(

        stats::prcomp(
          complete_data,
          center = TRUE,
          scale. = TRUE
        ),

        error = function(e) {

          if (verbose) {
            warning(
              "PCA failed: ",
              conditionMessage(e)
            )
          }

          NULL
        }
      )


      if (is.null(pca_model)) {
        pca_failed <- TRUE
      }
    }


    if (!pca_failed) {

      pca_scores <- as.data.frame(
        pca_model$x,
        check.names = FALSE
      )


      available_axes <- intersect(
        pca_axes,
        colnames(pca_scores)
      )


      unavailable_axes <- setdiff(
        pca_axes,
        available_axes
      )


      if (
        length(unavailable_axes) > 0 &&
        verbose
      ) {
        warning(
          "Requested PCA axes were not available: ",
          paste(
            unavailable_axes,
            collapse = ", "
          )
        )
      }


      if (length(available_axes) > 0) {

        pca_results <- pca_scores[
          ,
          available_axes,
          drop = FALSE
        ]


        for (axis in available_axes) {

          analysis_data[[axis]] <- NA_real_

          analysis_data[
            rownames(pca_results),
            axis
          ] <- pca_results[[axis]]
        }
      }
    }
  }


  # ---------------------------------------------------------------------------
  # 3. Helper functions
  # ---------------------------------------------------------------------------

  get_significance <- function(p) {

    if (is.na(p)) {
      return("")
    }

    if (p < sig_levels[1]) {

      "***"

    } else if (p < sig_levels[2]) {

      "**"

    } else if (p < sig_levels[3]) {

      "*"

    } else {

      "ns"
    }
  }


  calculate_signal <- function(tree,
                               trait,
                               method) {

    result <- tryCatch(

      {

        if (method == "lambda") {

          phytools::phylosig(
            tree,
            trait,
            method = "lambda",
            test = TRUE
          )

        } else {

          phytools::phylosig(
            tree,
            trait,
            method = "K",
            test = TRUE,
            nsim = nsim
          )
        }
      },

      error = function(e) {

        if (verbose) {
          warning(
            "Phylogenetic signal calculation failed for ",
            method,
            ": ",
            conditionMessage(e)
          )
        }

        NULL
      }
    )


    if (is.null(result)) {
      return(
        list(
          signal = NA_real_,
          p = NA_real_
        )
      )
    }


    signal <- if (method == "lambda") {
      result$lambda
    } else {
      result$K
    }


    list(
      signal = as.numeric(signal),
      p = as.numeric(result$P)
    )
  }


  # ---------------------------------------------------------------------------
  # 4. Calculate phylogenetic signal within each community
  # ---------------------------------------------------------------------------

  trait_names <- colnames(
    analysis_data
  )


  result_list <- vector(
    "list",
    nrow(com) *
      length(trait_names) *
      length(methods)
  )


  result_i <- 0L

  small_pool_count <- 0L


  for (i in seq_len(nrow(com))) {

    plot_name <- plot_names[i]


    plot_abundance <- unlist(
      com[i, , drop = TRUE],
      use.names = TRUE
    )


    present <- !is.na(plot_abundance) &
      plot_abundance > min_abundance


    present_species <- names(
      plot_abundance
    )[present]


    # Total number of species represented in the community
    n_sp_in_plot <- length(
      present_species
    )


    for (trait_name in trait_names) {

      # -----------------------------------------------------------------------
      # Trait values for species present in this community
      # -----------------------------------------------------------------------

      trait_values <- rep(
        NA_real_,
        n_sp_in_plot
      )

      names(trait_values) <- present_species


      trait_rows <- match(
        present_species,
        rownames(analysis_data)
      )


      matched <- !is.na(
        trait_rows
      )


      trait_values[matched] <- analysis_data[
        trait_rows[matched],
        trait_name
      ]


      # -----------------------------------------------------------------------
      # Trait coverage
      # -----------------------------------------------------------------------

      if (n_sp_in_plot > 0) {

        coverage <- paste0(
          round(
            mean(
              !is.na(trait_values)
            ) * 100,
            2
          ),
          " %"
        )

      } else {

        coverage <- NA_character_
      }


      # -----------------------------------------------------------------------
      # Species eligible for phylogenetic signal analysis
      # -----------------------------------------------------------------------

      analytical_species <- present_species[
        !is.na(trait_values) &
          present_species %in% phylo_tree$tip.label
      ]


      # Actual analytical sample size
      n_sp <- length(
        analytical_species
      )


      can_analyze <- n_sp >= 4


      if (!can_analyze) {
        small_pool_count <- small_pool_count + 1L
      }


      if (can_analyze) {

        final_tree <- ape::drop.tip(
          phylo_tree,
          setdiff(
            phylo_tree$tip.label,
            analytical_species
          )
        )


        final_trait <- trait_values[
          final_tree$tip.label
        ]


        if (
          length(
            unique(
              final_trait
            )
          ) < 2
        ) {

          can_analyze <- FALSE

          if (verbose) {
            warning(
              "Trait `",
              trait_name,
              "` has no variation in community `",
              plot_name,
              "` after matching with the phylogeny."
            )
          }
        }
      }


      # -----------------------------------------------------------------------
      # Calculate requested phylogenetic signal metrics
      # -----------------------------------------------------------------------

      for (method in methods) {

        result_i <- result_i + 1L


        if (can_analyze) {

          signal_result <- calculate_signal(
            tree = final_tree,
            trait = final_trait,
            method = method
          )

        } else {

          signal_result <- list(
            signal = NA_real_,
            p = NA_real_
          )
        }


        result_list[[result_i]] <- data.frame(
          plot = plot_name,
          trait = trait_name,
          coverage = coverage,
          n_sp = n_sp,
          signal = signal_result$signal,
          p = signal_result$p,
          significance = get_significance(
            signal_result$p
          ),
          method = method,
          n_sp_in_plot = n_sp_in_plot,
          stringsAsFactors = FALSE
        )
      }
    }
  }


  # ---------------------------------------------------------------------------
  # 5. Return results
  # ---------------------------------------------------------------------------

  results <- do.call(
    rbind,
    result_list[
      seq_len(result_i)
    ]
  )


  rownames(results) <- NULL


  analyzed_plots <- unique(
    results$plot[
      !is.na(results$signal)
    ]
  )


  if (
    small_pool_count > 0 &&
    verbose
  ) {
    warning(
      small_pool_count,
      " community-trait combination",
      if (small_pool_count == 1) "" else "s",
      " were not analyzed because fewer than four species had both ",
      "trait data and phylogenetic information."
    )
  }


  attr(
    results,
    "methods"
  ) <- methods

  attr(
    results,
    "pca_axes"
  ) <- pca_axes

  attr(
    results,
    "pca_failed"
  ) <- pca_failed

  attr(
    results,
    "sig_levels"
  ) <- sig_levels

  attr(
    results,
    "min_abundance"
  ) <- min_abundance

  attr(
    results,
    "nsim"
  ) <- nsim

  attr(
    results,
    "total_plots"
  ) <- nrow(com)

  attr(
    results,
    "analyzed_plots"
  ) <- length(
    analyzed_plots
  )


  if (!is.null(pca_results)) {
    attr(
      results,
      "pca_results"
    ) <- pca_results
  }


  if (!is.null(pca_model)) {
    attr(
      results,
      "pca_model"
    ) <- pca_model
  }


  results
}
