#' Estimate Phylogenetic Signal in Trait Data
#'
#' `pnc()` estimates phylogenetic signal in functional traits for a given
#' species pool. The function matches trait data with a phylogenetic tree,
#' optionally summarizes multivariate trait variation using principal
#' component analysis (PCA), and estimates phylogenetic signal using
#' Pagel's lambda or Blomberg's K.
#'
#' Phylogenetic signal describes the tendency for related species to resemble
#' one another in trait values. It can provide evidence relevant to
#' phylogenetic niche conservatism, but does not by itself constitute a
#' direct test of phylogenetic niche conservatism.
#'
#' @param trait_data A data frame or matrix containing continuous trait data,
#'   with species as rows and traits as columns. Species names must be stored
#'   as row names.
#' @param phylo_tree A phylogenetic tree object of class `"phylo"`.
#' @param methods Character vector specifying the phylogenetic signal metrics
#'   to calculate. Available options are `"lambda"` and `"K"`.
#' @param pca_axes Character vector specifying PCA axes to include, for example
#'   `c("PC1", "PC2")`. Set to `NULL` to skip PCA.
#' @param sig_levels Numeric vector of three significance thresholds used to
#'   assign significance symbols.
#' @param nsim Number of randomizations used for significance testing of
#'   Blomberg's K.
#' @param verbose Logical. If `TRUE`, warnings from PCA or phylogenetic signal
#'   estimation are reported.
#'
#' @return A data frame with one row for each trait-method combination:
#'   \itemize{
#'     \item `trait`: trait name or PCA axis.
#'     \item `coverage`: percentage of species with available trait values.
#'     \item `n_sp`: number of species actually included in the phylogenetic
#'       signal analysis after matching trait data with the phylogeny.
#'     \item `signal`: estimated phylogenetic signal.
#'     \item `p`: P value.
#'     \item `significance`: significance symbol.
#'     \item `method`: phylogenetic signal metric.
#'   }
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
#'   traits = c("LA", "LMA", "LeafN", "PlantHeight", "SeedMass", "SSD")
#' )
#'
#' # Pagel's lambda with PCA
#' pnc(
#'   subtraits,
#'   BCI$phy_species,
#'   methods = "lambda"
#' )
#'
#' # Pagel's lambda without PCA
#' pnc(
#'   subtraits,
#'   BCI$phy_species,
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

  trait_data <- as.data.frame(trait_data, check.names = FALSE)

  if (nrow(trait_data) == 0 || ncol(trait_data) == 0) {
    stop("`trait_data` must contain at least one species and one trait.")
  }

  if (anyDuplicated(rownames(trait_data))) {
    stop("Species names in `trait_data` must be unique.")
  }

  if (anyDuplicated(phylo_tree$tip.label)) {
    stop("Tip labels in `phylo_tree` must be unique.")
  }

  if (!inherits(phylo_tree, "phylo")) {
    stop("`phylo_tree` must be an object of class `phylo`.")
  }

  numeric_traits <- vapply(trait_data, is.numeric, logical(1))

  if (!all(numeric_traits)) {
    stop(
      "All trait columns must be numeric. Non-numeric traits: ",
      paste(names(numeric_traits)[!numeric_traits], collapse = ", ")
    )
  }

  if (length(intersect(rownames(trait_data), phylo_tree$tip.label)) == 0) {
    stop("No species are shared between `trait_data` and `phylo_tree`.")
  }

  methods <- match.arg(
    methods,
    choices = c("lambda", "K"),
    several.ok = TRUE
  )

  if (!is.null(pca_axes) && !is.character(pca_axes)) {
    stop("`pca_axes` must be a character vector or `NULL`.")
  }

  pca_axes <- unique(pca_axes)

  if (
    length(sig_levels) != 3 ||
    any(!is.finite(sig_levels)) ||
    any(diff(sig_levels) <= 0)
  ) {
    stop("`sig_levels` must contain three increasing significance thresholds.")
  }

  if (
    length(nsim) != 1 ||
    !is.finite(nsim) ||
    nsim < 1
  ) {
    stop("`nsim` must be a positive integer.")
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

    complete_rows <- stats::complete.cases(trait_data)
    complete_data <- trait_data[complete_rows, , drop = FALSE]

    if (nrow(complete_data) < 4) {

      pca_failed <- TRUE

      if (verbose) {
        warning(
          "PCA was skipped because fewer than four species ",
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

      if (length(unavailable_axes) > 0 && verbose) {
        warning(
          "Requested PCA axes were not available: ",
          paste(unavailable_axes, collapse = ", ")
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

    if (p <= sig_levels[1]) {
      "***"
    } else if (p <= sig_levels[2]) {
      "**"
    } else if (p <= sig_levels[3]) {
      "*"
    } else {
      "ns"
    }
  }


  calculate_signal <- function(tree, trait, method) {

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
      signal = signal,
      p = result$P
    )
  }


  # ---------------------------------------------------------------------------
  # 4. Calculate phylogenetic signal
  # ---------------------------------------------------------------------------

  trait_names <- colnames(analysis_data)

  result_list <- vector(
    "list",
    length(trait_names) * length(methods)
  )

  result_i <- 0L

  for (trait_name in trait_names) {

    trait_values <- analysis_data[[trait_name]]

    coverage <- paste0(
      round(
        mean(!is.na(trait_values)) * 100,
        2
      ),
      " %"
    )

    valid <- !is.na(trait_values)

    trait_vector <- trait_values[valid]
    species <- rownames(analysis_data)[valid]

    names(trait_vector) <- species

    common_species <- intersect(
      species,
      phylo_tree$tip.label
    )

    # Actual analytical sample size
    n_sp <- length(common_species)

    can_analyze <- n_sp >= 4

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

      if (length(unique(trait_for_analysis)) < 2) {

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
          pruned_tree,
          trait_for_analysis,
          method
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
        significance = get_significance(signal_result$p),
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

  attr(results, "methods") <- methods
  attr(results, "pca_axes") <- pca_axes
  attr(results, "pca_failed") <- pca_failed
  attr(results, "sig_levels") <- sig_levels
  attr(results, "nsim") <- nsim

  if (!is.null(pca_results)) {
    attr(results, "pca_results") <- pca_results
  }

  if (!is.null(pca_model)) {
    attr(results, "pca_model") <- pca_model
  }

  results
}
