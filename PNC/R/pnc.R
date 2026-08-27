#' Estimate Phylogenetic Signal in Quantitative Traits
#'
#' `pnc()` estimates phylogenetic signal in quantitative traits for a focal
#' set of taxa. The function matches trait data to a phylogenetic tree and
#' estimates Pagel's lambda, Blomberg's K, or both for individual traits and,
#' when requested, principal component (PC) axes.
#'
#' Phylogenetic signal describes the tendency for related taxa to resemble
#' one another in measured trait values. The resulting estimates can be used
#' to evaluate whether phylogenetic relatedness is informative about trait
#' similarity in the analyzed taxon pool. Their ecological interpretation
#' depends on the traits supplied and the taxa retained in the analysis.
#'
#' @param trait_data A data frame or matrix containing quantitative trait
#'   values, with taxa as rows and traits as columns. Taxon names must be
#'   stored as unique row names and are matched exactly to the tip labels of
#'   `phylo_tree`. Trait columns must be numeric. Missing values may be `NA`;
#'   all other values must be finite.
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
#' @param nsim Number of randomizations used to test Blomberg's K. This
#'   argument is ignored when only Pagel's lambda is requested. The default
#'   is 1000.
#' @param verbose Logical. If `TRUE`, warnings about skipped or failed PCA and
#'   phylogenetic signal analyses are reported.
#'
#' @return A data frame with one row for each trait-method combination:
#'   \itemize{
#'     \item `trait`: name of the individual trait or PC axis.
#'     \item `coverage`: formatted percentage of input taxa with an available
#'       value. For an individual trait, this is trait availability before
#'       phylogenetic matching. For a PC axis, it represents complete-case
#'       coverage among taxa matched to the phylogeny.
#'     \item `n_sp`: number of taxa included after matching non-missing trait
#'       values to the phylogeny.
#'     \item `signal`: estimated Pagel's lambda or Blomberg's K.
#'     \item `p`: P-value for the phylogenetic signal estimate.
#'     \item `significance`: significance symbol assigned using `sig_levels`.
#'     \item `method`: phylogenetic signal metric used.
#'   }
#'
#'   The returned data frame also stores `methods`, `pca_axes`, `pca_failed`,
#'   `sig_levels`, and `nsim` as attributes. When PCA is successfully fitted,
#'   the PCA scores and fitted model are stored in the `pca_results` and
#'   `pca_model` attributes, respectively.
#'
#' @details
#' Individual traits are analyzed separately. Missing values in one trait
#' therefore do not alter the taxa included in the analysis of another trait.
#'
#' When PC axes are requested, PCA is performed once using taxa represented
#' in both `trait_data` and `phylo_tree` and having complete observations
#' across all supplied traits. Trait values are centered and scaled before
#' PCA. The requested PC scores are then analyzed in the same way as
#' individual traits.
#'
#' PCA is skipped if fewer than four matched taxa have complete trait data or
#' if the PCA cannot be fitted. This does not prevent the individual traits
#' from being analyzed.
#'
#' For each individual trait or PC axis, phylogenetic signal is estimated only
#' when at least four taxa have both an observed value and a matching
#' phylogenetic tip, and when the retained values show variation. Otherwise,
#' `signal` and `p` are returned as `NA`. The four-taxon threshold is a
#' computational minimum used by the function, not a general recommendation
#' for adequate sampling.
#'
#' @examples
#' \donttest{
#' data(BCI)
#' data(TRY)
#'
#' sp <- colnames(BCI$com)
#'
#' subtraits <- extract_traits(
#'   sp,
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
#' # Pagel's lambda for individual traits and PC1 and PC2
#' pnc(
#'   subtraits,
#'   BCI$phy_species,
#'   methods = "lambda"
#' )
#'
#' # Pagel's lambda for individual traits only
#' pnc(
#'   subtraits,
#'   BCI$phy_species,
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
pnc <- function(trait_data,
                phylo_tree,
                methods = "lambda",
                pca_axes = c("PC1", "PC2"),
                sig_levels = c(0.001, 0.01, 0.05),
                nsim = 1000,
                verbose = TRUE) {

  # ---------------------------------------------------------------------------
  # 1. Check inputs
  # ---------------------------------------------------------------------------

  if (!is.data.frame(trait_data) && !is.matrix(trait_data)) {
    stop("`trait_data` must be a data frame or matrix.")
  }

  if (!inherits(phylo_tree, "phylo")) {
    stop("`phylo_tree` must be an object of class `phylo`.")
  }

  trait_data <- as.data.frame(
    trait_data,
    check.names = FALSE
  )

  if (nrow(trait_data) == 0 || ncol(trait_data) == 0) {
    stop("`trait_data` must contain at least one species and one trait.")
  }

  if (is.null(rownames(trait_data))) {
    stop("`trait_data` must have species names as row names.")
  }

  if (anyDuplicated(rownames(trait_data))) {
    stop("Species names in `trait_data` must be unique.")
  }

  if (anyDuplicated(phylo_tree$tip.label)) {
    stop("Tip labels in `phylo_tree` must be unique.")
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

  reference_species <- intersect(
    rownames(trait_data),
    phylo_tree$tip.label
  )

  if (length(reference_species) == 0) {
    stop(
      "No species are shared between `trait_data` and `phylo_tree`."
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

    pca_axes <- unique(pca_axes)

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
    stop("`nsim` must be a positive integer.")
  }

  if (
    length(verbose) != 1 ||
    !is.logical(verbose) ||
    is.na(verbose)
  ) {
    stop("`verbose` must be TRUE or FALSE.")
  }

  nsim <- as.integer(nsim)


  # ---------------------------------------------------------------------------
  # 2. Add PCA scores
  # ---------------------------------------------------------------------------

  analysis_data <- trait_data

  pca_model <- NULL
  pca_results <- NULL
  pca_failed <- FALSE

  if (length(pca_axes) > 0) {

    # PCA is based only on species represented in both
    # the trait data and the phylogeny.
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
          "shared by `trait_data` and `phylo_tree` had complete trait data."
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
  # 4. Calculate phylogenetic signal
  # ---------------------------------------------------------------------------

  trait_names <- colnames(
    analysis_data
  )

  result_list <- vector(
    "list",
    length(trait_names) *
      length(methods)
  )

  result_i <- 0L

  for (trait_name in trait_names) {

    trait_values <- analysis_data[[trait_name]]

    coverage <- paste0(
      round(
        mean(
          !is.na(trait_values)
        ) * 100,
        2
      ),
      " %"
    )

    valid <- !is.na(
      trait_values
    )

    trait_vector <- trait_values[
      valid
    ]

    species <- rownames(
      analysis_data
    )[valid]

    names(trait_vector) <- species

    common_species <- intersect(
      species,
      phylo_tree$tip.label
    )

    # Number of species actually entering
    # the phylogenetic signal analysis.
    n_sp <- length(
      common_species
    )

    can_analyze <- n_sp >= 4

    if (!can_analyze && verbose) {

      warning(
        "Trait `",
        trait_name,
        "` was not analyzed because fewer than four species ",
        "had both trait data and phylogenetic information."
      )
    }

    if (can_analyze) {

      pruned_tree <- ape::drop.tip(
        phylo_tree,
        setdiff(
          phylo_tree$tip.label,
          common_species
        )
      )

      trait_for_analysis <- trait_vector[
        pruned_tree$tip.label
      ]

      if (
        length(
          unique(
            trait_for_analysis
          )
        ) < 2
      ) {

        can_analyze <- FALSE

        if (verbose) {
          warning(
            "Trait `",
            trait_name,
            "` has no variation after matching with the phylogeny."
          )
        }
      }
    }

    for (method in methods) {

      result_i <- result_i + 1L

      if (can_analyze) {

        signal_result <- calculate_signal(
          tree = pruned_tree,
          trait = trait_for_analysis,
          method = method
        )

      } else {

        signal_result <- list(
          signal = NA_real_,
          p = NA_real_
        )
      }

      result_list[[result_i]] <- data.frame(
        trait = trait_name,
        coverage = coverage,
        n_sp = n_sp,
        signal = signal_result$signal,
        p = signal_result$p,
        significance = get_significance(
          signal_result$p
        ),
        method = method,
        stringsAsFactors = FALSE
      )
    }
  }


  # ---------------------------------------------------------------------------
  # 5. Return results
  # ---------------------------------------------------------------------------

  results <- do.call(
    rbind,
    result_list
  )

  rownames(results) <- NULL

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
    "nsim"
  ) <- nsim

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
