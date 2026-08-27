#' Estimate Phylogenetic Signal Across Multiple Communities
#'
#' `compnc()` estimates phylogenetic signal in quantitative traits across
#' multiple ecological communities. The function identifies the species
#' present in each community, matches them to a trait matrix and phylogenetic
#' tree, and estimates Pagel's lambda, Blomberg's K, or both within each
#' resulting local species pool.
#'
#' When PC axes are requested, one principal component analysis (PCA) is
#' fitted using the pooled set of eligible species across all communities.
#' The resulting species scores are then used in every community, so that a
#' given PC axis represents the same dimension of trait variation across
#' communities.
#'
#' Phylogenetic signal describes the tendency for related species to resemble
#' one another in measured trait values. The community-specific estimates can
#' be used to evaluate whether phylogenetic relatedness is informative about
#' trait similarity within each local species pool and to compare this
#' relationship among communities.
#'
#' @param com A numeric community matrix or data frame with communities as
#'   rows and species as columns. Cell values represent species abundance or
#'   occurrence and must be non-negative. Species names must be stored as
#'   unique column names. Row names are used as community names; when they are
#'   absent, names are generated automatically. Missing values are treated as
#'   absence.
#' @param trait_data A data frame or matrix containing quantitative trait
#'   values, with species as rows and traits as columns. Species names must be
#'   stored as unique row names and are matched exactly to the species names
#'   in `com` and the tip labels of `phylo_tree`. Trait columns must be
#'   numeric. Missing values may be `NA`; all other values must be finite.
#' @param phylo_tree A phylogenetic tree object of class `"phylo"` with unique
#'   tip labels.
#' @param methods Character vector specifying the phylogenetic signal metrics
#'   to calculate. Available options are `"lambda"` and `"K"`. The default is
#'   `"lambda"`.
#' @param pca_axes Character vector naming the PC axes to include in the
#'   phylogenetic signal analysis. The default is `c("PC1", "PC2")`. Set to
#'   `NULL` to skip PCA.
#' @param sig_levels Numeric vector of three increasing P-value thresholds
#'   used to assign significance symbols. The default is
#'   `c(0.001, 0.01, 0.05)`.
#' @param min_abundance Minimum abundance required for a species to be treated
#'   as present in a community. A species is considered present when its value
#'   is greater than `min_abundance`. The default is 0.
#' @param nsim Number of randomizations used to test Blomberg's K. This
#'   argument is ignored when only Pagel's lambda is requested. The default
#'   is 1000.
#' @param verbose Logical. If `TRUE`, warnings about skipped or failed PCA and
#'   phylogenetic signal analyses are reported.
#'
#' @return A data frame with one row for each
#'   community-trait-method combination:
#'   \itemize{
#'     \item `plot`: community name.
#'     \item `trait`: name of the individual trait or PC axis.
#'     \item `coverage`: formatted percentage of species present in the
#'       community with an available value. For an individual trait, this is
#'       trait availability before phylogenetic matching. For a PC axis, it is
#'       the percentage of species with a score from the pooled PCA.
#'     \item `n_sp`: number of species included after matching the trait or PC
#'       score to the phylogeny.
#'     \item `signal`: estimated Pagel's lambda or Blomberg's K.
#'     \item `p`: P-value for the phylogenetic signal estimate.
#'     \item `significance`: significance symbol assigned using `sig_levels`.
#'     \item `method`: phylogenetic signal metric used.
#'     \item `n_sp_in_plot`: number of species present in the original local
#'       community according to `min_abundance`, before trait and phylogenetic
#'       matching.
#'   }
#'
#'   The returned data frame also stores `methods`, `pca_axes`, `pca_failed`,
#'   `sig_levels`, `min_abundance`, `nsim`, `total_plots`, and
#'   `analyzed_plots` as attributes. When PCA is successfully fitted, the
#'   pooled species scores and fitted model are stored in the `pca_results`
#'   and `pca_model` attributes, respectively.
#'
#' @details
#' Individual traits are analyzed separately. Missing values in one trait
#' therefore do not alter the species included in the analysis of another
#' trait.
#'
#' When PC axes are requested, the pooled PCA includes species that occur in
#' at least one community, are represented in both `trait_data` and
#' `phylo_tree`, and have complete observations across all supplied traits.
#' Trait values are centered and scaled before PCA. The PCA is fitted once,
#' and the same species scores are used in every community.
#'
#' The output reports both `n_sp_in_plot` and trait-specific `n_sp`, allowing
#' users to compare the original local pool with the species retained after
#' trait and phylogenetic matching.
#'
#' Phylogenetic signal is estimated only when at least four matched species
#' have an observed value and the retained values show variation. Otherwise,
#' `signal` and `p` are returned as `NA`. The four-species threshold is a
#' computational minimum used by the function, not a general recommendation
#' for adequate sampling.
#'
#' Community-specific estimates may be sensitive to the number and
#' phylogenetic composition of the retained species. Results from small or
#' phylogenetically restricted local pools should therefore be interpreted
#' together with the reported `n_sp` and the species retained in each
#' analysis.
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
#' # Pagel's lambda using one pooled PCA across all communities
#' compnc(
#'   HimalayanBirds$com,
#'   subtraits,
#'   HimalayanBirds$phy_species,
#'   methods = "lambda"
#' )
#'
#' # Pagel's lambda for individual traits only
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
#' Pagel, M. (1999). Inferring the historical patterns of biological
#' evolution. Nature, 401, 877-884.
#' \doi{10.1038/44766}
#'
#' Blomberg, S. P., Garland, T., Jr., and Ives, A. R. (2003). Testing for
#' phylogenetic signal in comparative data: Behavioral traits are more
#' labile. Evolution, 57, 717-745.
#' \doi{10.1111/j.0014-3820.2003.tb00285.x}
#'
#' Münkemüller, T., Lavergne, S., Bzeznik, B., Dray, S., Jombart, T.,
#' Schiffers, K., and Thuiller, W. (2012). How to measure and test
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
