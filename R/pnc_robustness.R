#' Assess Sensitivity of Phylogenetic Signal to Missing Trait Data
#'
#' `pnc_robustness()` evaluates how the observed pattern of missing trait data
#' affects inference of phylogenetic signal based on Pagel's lambda.
#'
#' For each trait, the observed lambda is first estimated using `pnc()`.
#' Complete trait data are then simulated under Brownian motion on a reference
#' phylogeny transformed according to the observed lambda. Each simulated trait
#' is analyzed twice: first using all species in the reference species pool,
#' and then after applying the observed missing-data pattern of the focal trait.
#'
#' The paired comparison isolates the effect of trait-data incompleteness on
#' phylogenetic signal estimation.
#'
#' @param trait_data A data frame or matrix containing continuous trait data,
#'   with species as rows and traits as columns. Species names must be stored
#'   as row names.
#' @param phylo_tree A phylogenetic tree object of class `"phylo"` with branch
#'   lengths.
#' @param n_simulations Number of paired simulations performed for each trait.
#'   Default is 100.
#' @param alpha_level Significance threshold used to compare statistical
#'   conclusions between complete and incomplete simulated data.
#'   Default is 0.05.
#' @param verbose Logical. If `TRUE`, a progress bar and warnings are shown.
#'   Default is `TRUE`.
#'
#' @return A data frame containing the observed phylogenetic signal results
#'   from `pnc()` and the following sensitivity metrics:
#'   \itemize{
#'     \item `simulation_lambda`: lambda value used to generate simulated trait
#'       data. This equals the observed lambda unless the latter exceeds the
#'       maximum supported by the reference phylogeny.
#'     \item `consistency`: percentage of successful paired simulations in which
#'       complete and incomplete data produced the same significance conclusion.
#'     \item `signal_bias`: mean paired difference between lambda estimated
#'       from incomplete and complete simulated data
#'       (`lambda_missing - lambda_complete`).
#'     \item `signal_sd`: standard deviation of the paired lambda differences.
#'     \item `n_successful`: number of simulations in which both complete and
#'       incomplete analyses were successfully estimated.
#'   }
#'
#' @details
#' The analysis is performed separately for each original trait. The reference
#' species pool consists of species represented in both `trait_data` and
#' `phylo_tree`. Species without observed values for the focal trait are retained
#' in the complete simulated data and removed only when the empirical
#' missing-data pattern is applied.
#'
#' For each trait, simulation is conditional on the Pagel's lambda estimated
#' from the observed data. The reference phylogeny is transformed according to
#' this lambda, and trait values are simulated under Brownian motion on the
#' transformed tree.
#'
#' The maximum valid lambda depends on the branch-length structure of the
#' reference phylogeny. An observed lambda estimated from a trait-specific
#' pruned tree may occasionally exceed the maximum value that can be applied
#' to the larger reference phylogeny. In this case, the largest valid lambda
#' for the reference phylogeny is used only for simulation. The observed lambda
#' estimate itself is not modified.
#'
#' Each simulated trait is analyzed before and after applying the observed
#' missing-data pattern. Because both analyses use the same simulated trait,
#' their paired difference reflects the effect of trait-data incompleteness
#' rather than stochastic differences among simulations.
#'
#' Simulated lambda estimates are not forced to equal the observed lambda.
#' Their natural stochastic variation is retained.
#'
#' Repeated Pagel's lambda estimation is accelerated by caching the
#' phylogeny-specific covariance structure and using an algebraically
#' equivalent spectral representation of the likelihood. If the fast
#' calculation fails, `phytools::phylosig()` is used as a fallback.
#'
#' PCA axes are not included because evaluating their missing-data sensitivity
#' would require multivariate trait simulation that preserves covariance among
#' traits together with trait-specific missing-data patterns.
#'
#' Blomberg's K is not used as the generating parameter because K is a
#' Brownian-motion-referenced summary statistic rather than a direct parameter
#' of the trait-generating model.
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
#'   traits = c("LA", "LMA", "LeafN", "PlantHeight")
#' )
#'
#' set.seed(123)
#'
#' pnc_robustness(
#'   subtraits,
#'   BCI$phy_species,
#'   n_simulations = 10
#' )
#' }
#'
#' @references
#' Pagel, M. (1999). Inferring the historical patterns of biological
#' evolution. Nature, 401, 877-884.
#'
#' Münkemüller, T., Lavergne, S., Bzeznik, B., Dray, S., Jombart, T.,
#' Schiffers, K. and Thuiller, W. (2012). How to measure and test
#' phylogenetic signal. Methods in Ecology and Evolution, 3, 743-756.
#' \doi{10.1111/j.2041-210X.2012.00196.x}
#'
#' @export
#' @import geiger
pnc_robustness <- function(trait_data,
                           phylo_tree,
                           n_simulations = 100,
                           alpha_level = 0.05,
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

  trait_species <- rownames(trait_data)

  if (
    is.null(trait_species) ||
    identical(
      trait_species,
      as.character(seq_len(nrow(trait_data)))
    )
  ) {
    stop("`trait_data` must have species names as row names.")
  }

  trait_data <- as.data.frame(
    trait_data,
    check.names = FALSE
  )

  if (nrow(trait_data) == 0 || ncol(trait_data) == 0) {
    stop(
      "`trait_data` must contain at least one species and one trait."
    )
  }

  if (anyDuplicated(rownames(trait_data))) {
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

  if (is.null(phylo_tree$edge.length)) {
    stop(
      "`phylo_tree` must contain branch lengths."
    )
  }

  if (
    any(!is.finite(phylo_tree$edge.length)) ||
    any(phylo_tree$edge.length < 0)
  ) {
    stop(
      "`phylo_tree` contains invalid branch lengths."
    )
  }

  if (
    length(n_simulations) != 1 ||
    !is.finite(n_simulations) ||
    n_simulations < 1 ||
    n_simulations %% 1 != 0
  ) {
    stop(
      "`n_simulations` must be a positive integer."
    )
  }

  if (
    length(alpha_level) != 1 ||
    !is.finite(alpha_level) ||
    alpha_level <= 0 ||
    alpha_level >= 1
  ) {
    stop(
      "`alpha_level` must be between 0 and 1."
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

  n_simulations <- as.integer(
    n_simulations
  )


  # ---------------------------------------------------------------------------
  # 2. Define reference species pool
  # ---------------------------------------------------------------------------

  reference_species <- intersect(
    rownames(trait_data),
    phylo_tree$tip.label
  )

  if (length(reference_species) < 4) {
    stop(
      "Fewer than four species are shared between `trait_data` ",
      "and `phylo_tree`."
    )
  }

  reference_tree <- ape::drop.tip(
    phylo_tree,
    setdiff(
      phylo_tree$tip.label,
      reference_species
    )
  )

  # Use phylogenetic tip order throughout the simulation.
  reference_species <- reference_tree$tip.label


  # ---------------------------------------------------------------------------
  # 3. Calculate observed phylogenetic signal
  # ---------------------------------------------------------------------------

  observed_results <- pnc(
    trait_data = trait_data,
    phylo_tree = phylo_tree,
    methods = "lambda",
    pca_axes = NULL,
    verbose = FALSE
  )


  # ---------------------------------------------------------------------------
  # 4. Helper functions
  # ---------------------------------------------------------------------------

  # Standard Pagel's lambda estimator used as a fallback.
  estimate_lambda_standard <- function(tree,
                                       trait) {

    trait <- trait[
      tree$tip.label
    ]

    result <- tryCatch(
      phytools::phylosig(
        tree,
        trait,
        method = "lambda",
        test = TRUE
      ),
      error = function(e) NULL
    )

    if (
      is.null(result) ||
      !is.finite(result$lambda) ||
      !is.finite(result$P)
    ) {
      return(
        list(
          lambda = NA_real_,
          p = NA_real_
        )
      )
    }

    list(
      lambda = as.numeric(result$lambda),
      p = as.numeric(result$P)
    )
  }


  # Prepare phylogeny-specific quantities for repeated lambda estimation.
  prepare_lambda_fast <- function(tree) {

    C <- ape::vcv.phylo(
      tree
    )

    d <- diag(C)

    if (
      any(!is.finite(d)) ||
      any(d <= 0)
    ) {
      stop(
        "Invalid diagonal values in the phylogenetic covariance matrix."
      )
    }

    inv_sqrt_d <- 1 / sqrt(d)

    R <- C * outer(
      inv_sqrt_d,
      inv_sqrt_d
    )

    eig <- eigen(
      R,
      symmetric = TRUE
    )

    u <- as.numeric(
      crossprod(
        eig$vectors,
        inv_sqrt_d
      )
    )

    max_lambda <- utils::getFromNamespace(
      "maxLambda",
      "phytools"
    )(tree)

    if (
      length(max_lambda) != 1 ||
      !is.finite(max_lambda) ||
      max_lambda <= 0
    ) {
      stop(
        "Could not determine a valid maximum Pagel's lambda."
      )
    }

    list(
      tip.label = rownames(C),
      eigenvalues = eig$values,
      eigenvectors = eig$vectors,
      inv_sqrt_d = inv_sqrt_d,
      u = u,
      logdet_D = sum(log(d)),
      max_lambda = max_lambda
    )
  }


  # Fast Pagel's lambda estimator using cached spectral decomposition.
  fit_lambda_fast <- function(prepared,
                              trait,
                              niter = 10L) {

    result <- tryCatch(
      {

        trait <- trait[
          prepared$tip.label
        ]

        if (
          any(!is.finite(trait)) ||
          length(unique(trait)) < 2
        ) {
          stop(
            "Trait values are invalid."
          )
        }

        n <- length(
          trait
        )

        v <- as.numeric(
          crossprod(
            prepared$eigenvectors,
            trait * prepared$inv_sqrt_d
          )
        )

        eigenvalues <- prepared$eigenvalues
        u <- prepared$u
        logdet_D <- prepared$logdet_D


        likelihood_lambda <- function(lambda) {

          lambda_eigenvalues <-
            1 + lambda * (eigenvalues - 1)

          if (
            any(!is.finite(lambda_eigenvalues)) ||
            any(lambda_eigenvalues <= 0)
          ) {
            return(-Inf)
          }

          inverse_eigenvalues <-
            1 / lambda_eigenvalues

          denominator <- sum(
            u^2 *
              inverse_eigenvalues
          )

          if (
            !is.finite(denominator) ||
            denominator <= 0
          ) {
            return(-Inf)
          }

          ancestral_mean <-
            sum(
              u *
                v *
                inverse_eigenvalues
            ) /
            denominator

          residual_projection <-
            v - ancestral_mean * u

          quadratic_form <- sum(
            residual_projection^2 *
              inverse_eigenvalues
          )

          if (
            !is.finite(quadratic_form) ||
            quadratic_form <= 0
          ) {
            return(-Inf)
          }

          sigma2 <-
            quadratic_form / n

          if (
            !is.finite(sigma2) ||
            sigma2 <= 0
          ) {
            return(-Inf)
          }

          logdet_C_lambda <-
            logdet_D +
            sum(
              log(
                lambda_eigenvalues
              )
            )

          log_likelihood <-
            -quadratic_form / (2 * sigma2) -
            n * log(2 * pi) / 2 -
            (
              n * log(sigma2) +
                logdet_C_lambda
            ) / 2

          as.numeric(
            log_likelihood
          )
        }


        # Match the interval-search strategy used by phylosig().
        max_lambda <- prepared$max_lambda

        intervals <- cbind(
          seq(
            0,
            max_lambda - max_lambda / niter,
            length.out = niter
          ),
          seq(
            max_lambda / niter,
            max_lambda,
            length.out = niter
          )
        )

        fits <- lapply(
          seq_len(niter),
          function(i) {
            stats::optimize(
              f = likelihood_lambda,
              interval = intervals[i, ],
              maximum = TRUE
            )
          }
        )

        likelihoods <- vapply(
          fits,
          function(x) {
            x$objective
          },
          numeric(1)
        )

        best_fit <- fits[[
          which.max(
            likelihoods
          )
        ]]

        lambda_hat <-
          best_fit$maximum

        logL <-
          best_fit$objective

        logL0 <-
          likelihood_lambda(0)

        likelihood_ratio <-
          2 * (logL - logL0)

        p_value <- stats::pchisq(
          likelihood_ratio,
          df = 1,
          lower.tail = FALSE
        )

        list(
          lambda = lambda_hat,
          p = p_value
        )
      },

      error = function(e) NULL
    )

    if (
      is.null(result) ||
      !is.finite(result$lambda) ||
      !is.finite(result$p)
    ) {
      return(
        list(
          lambda = NA_real_,
          p = NA_real_
        )
      )
    }

    result
  }


  # Use the fast estimator and fall back to phylosig() if necessary.
  estimate_lambda <- function(tree,
                              prepared,
                              trait) {

    if (!is.null(prepared)) {

      fast_result <- fit_lambda_fast(
        prepared = prepared,
        trait = trait
      )

      if (
        is.finite(fast_result$lambda) &&
        is.finite(fast_result$p)
      ) {
        return(
          fast_result
        )
      }
    }

    estimate_lambda_standard(
      tree = tree,
      trait = trait
    )
  }


  # Apply Pagel's lambda transformation.

  rescale_lambda <- function(tree,
                             lambda) {

    transformed_tree <- tryCatch(

      withCallingHandlers(

        phytools::rescale(
          tree,
          model = "lambda",
          lambda = lambda
        ),

        warning = function(w) {

          if (
            grepl(
              "negative branch lengths encountered",
              conditionMessage(w),
              fixed = TRUE
            )
          ) {
            invokeRestart(
              "muffleWarning"
            )
          }
        }
      ),

      error = function(e) NULL
    )

    if (
      is.null(transformed_tree) ||
      !inherits(transformed_tree, "phylo")
    ) {
      return(NULL)
    }

    if (
      is.null(transformed_tree$edge.length) ||
      any(!is.finite(transformed_tree$edge.length)) ||
      any(transformed_tree$edge.length < 0)
    ) {
      return(NULL)
    }

    transformed_tree
  }

  # Find the largest lambda supported by the reference phylogeny.
  find_valid_lambda <- function(tree,
                                lambda_upper,
                                n_iterations = 60L) {

    lambda_lower <- 1

    tree_lower <- rescale_lambda(
      tree = tree,
      lambda = lambda_lower
    )

    if (is.null(tree_lower)) {
      return(
        list(
          lambda = NA_real_,
          tree = NULL
        )
      )
    }

    for (iteration in seq_len(n_iterations)) {

      lambda_mid <-
        (
          lambda_lower +
            lambda_upper
        ) / 2

      tree_mid <- rescale_lambda(
        tree = tree,
        lambda = lambda_mid
      )

      if (is.null(tree_mid)) {

        lambda_upper <-
          lambda_mid

      } else {

        lambda_lower <-
          lambda_mid

        tree_lower <-
          tree_mid
      }
    }

    list(
      lambda = lambda_lower,
      tree = tree_lower
    )
  }


  # Standard empty sensitivity result.
  empty_sensitivity <- function() {

    data.frame(
      simulation_lambda = NA_real_,
      consistency = NA_character_,
      signal_bias = NA_real_,
      signal_sd = NA_real_,
      n_successful = 0L
    )
  }


  # ---------------------------------------------------------------------------
  # 5. Prepare fast engine for reference tree
  # ---------------------------------------------------------------------------

  prepared_reference <- tryCatch(

    prepare_lambda_fast(
      reference_tree
    ),

    error = function(e) {

      if (verbose) {
        warning(
          "Fast Pagel's lambda preparation failed for the reference ",
          "phylogeny. Standard `phytools::phylosig()` will be used.",
          call. = FALSE
        )
      }

      NULL
    }
  )


  # ---------------------------------------------------------------------------
  # 6. Set up progress bar
  # ---------------------------------------------------------------------------

  total_steps <-
    nrow(observed_results) *
    n_simulations

  progress_step <- 0L

  pb <- NULL

  if (verbose) {

    pb <- utils::txtProgressBar(
      min = 0,
      max = total_steps,
      style = 3
    )

    on.exit(
      {
        if (!is.null(pb)) {
          close(pb)
        }
      },
      add = TRUE
    )
  }


  # ---------------------------------------------------------------------------
  # 7. Run paired simulations
  # ---------------------------------------------------------------------------

  sensitivity_results <- vector(
    "list",
    nrow(observed_results)
  )


  for (i in seq_len(nrow(observed_results))) {

    trait_name <-
      observed_results$trait[i]

    lambda_obs <-
      observed_results$signal[i]


    # -------------------------------------------------------------------------
    # Skip traits for which observed lambda could not be estimated
    # -------------------------------------------------------------------------

    if (!is.finite(lambda_obs)) {

      sensitivity_results[[i]] <-
        empty_sensitivity()

      progress_step <-
        progress_step +
        n_simulations

      if (verbose) {
        utils::setTxtProgressBar(
          pb,
          progress_step
        )
      }

      next
    }


    # -------------------------------------------------------------------------
    # Identify the empirical missing-data pattern
    # -------------------------------------------------------------------------

    trait_values <-
      trait_data[[trait_name]]

    names(trait_values) <-
      rownames(trait_data)


    observed_species <- reference_species[
      !is.na(
        trait_values[
          reference_species
        ]
      )
    ]


    if (length(observed_species) < 4) {

      sensitivity_results[[i]] <-
        empty_sensitivity()

      progress_step <-
        progress_step +
        n_simulations

      if (verbose) {
        utils::setTxtProgressBar(
          pb,
          progress_step
        )
      }

      next
    }


    observed_tree <- ape::drop.tip(
      reference_tree,
      setdiff(
        reference_species,
        observed_species
      )
    )


    # -------------------------------------------------------------------------
    # Prepare fast engine for incomplete-data analysis
    # -------------------------------------------------------------------------

    prepared_observed <- tryCatch(

      prepare_lambda_fast(
        observed_tree
      ),

      error = function(e) {

        if (verbose) {
          warning(
            "Fast Pagel's lambda preparation failed for trait `",
            trait_name,
            "`. Standard `phytools::phylosig()` will be used ",
            "for the incomplete-data analysis.",
            call. = FALSE
          )
        }

        NULL
      }
    )


    # -------------------------------------------------------------------------
    # Construct lambda-transformed reference phylogeny
    # -------------------------------------------------------------------------

    simulation_lambda <-
      lambda_obs

    simulation_tree <- rescale_lambda(
      tree = reference_tree,
      lambda = simulation_lambda
    )


    # An observed lambda > 1 estimated from a smaller trait-specific tree may
    # exceed the maximum supported by the larger reference phylogeny.
    if (
      is.null(simulation_tree) &&
      lambda_obs > 1
    ) {

      lambda_limit <- find_valid_lambda(
        tree = reference_tree,
        lambda_upper = lambda_obs
      )

      simulation_lambda <-
        lambda_limit$lambda

      simulation_tree <-
        lambda_limit$tree


      if (
        verbose &&
        is.finite(simulation_lambda)
      ) {

        warning(
          "Observed lambda for trait `",
          trait_name,
          "` (",
          formatC(
            lambda_obs,
            format = "f",
            digits = 6
          ),
          ") exceeds the maximum supported by the reference phylogeny. ",
          "Lambda = ",
          formatC(
            simulation_lambda,
            format = "f",
            digits = 6
          ),
          " was used for simulation.",
          call. = FALSE
        )
      }
    }


    # -------------------------------------------------------------------------
    # Skip if a valid simulation tree cannot be constructed
    # -------------------------------------------------------------------------

    if (is.null(simulation_tree)) {

      sensitivity_results[[i]] <-
        empty_sensitivity()

      progress_step <-
        progress_step +
        n_simulations

      if (verbose) {

        warning(
          "Could not construct a valid lambda-transformed reference ",
          "phylogeny for trait `",
          trait_name,
          "`.",
          call. = FALSE
        )

        utils::setTxtProgressBar(
          pb,
          progress_step
        )
      }

      next
    }


    # -------------------------------------------------------------------------
    # Storage for paired simulations
    # -------------------------------------------------------------------------

    lambda_difference <- rep(
      NA_real_,
      n_simulations
    )

    same_conclusion <- rep(
      NA,
      n_simulations
    )


    # -------------------------------------------------------------------------
    # Simulation loop
    # -------------------------------------------------------------------------

    for (sim_i in seq_len(n_simulations)) {

      simulated_trait <- tryCatch(

        phytools::fastBM(
          simulation_tree
        ),

        error = function(e) NULL
      )


      if (!is.null(simulated_trait)) {

        simulated_trait <- simulated_trait[
          reference_tree$tip.label
        ]


        # Complete-data analysis
        complete_result <- estimate_lambda(
          tree = reference_tree,
          prepared = prepared_reference,
          trait = simulated_trait
        )


        # Apply empirical missing-data pattern
        missing_trait <- simulated_trait[
          observed_tree$tip.label
        ]


        # Incomplete-data analysis
        missing_result <- estimate_lambda(
          tree = observed_tree,
          prepared = prepared_observed,
          trait = missing_trait
        )


        # Keep only successful paired analyses
        if (
          is.finite(complete_result$lambda) &&
          is.finite(complete_result$p) &&
          is.finite(missing_result$lambda) &&
          is.finite(missing_result$p)
        ) {

          lambda_difference[sim_i] <-
            missing_result$lambda -
            complete_result$lambda


          complete_significant <-
            complete_result$p <
            alpha_level

          missing_significant <-
            missing_result$p <
            alpha_level


          same_conclusion[sim_i] <-
            complete_significant ==
            missing_significant
        }
      }


      progress_step <-
        progress_step +
        1L

      if (verbose) {
        utils::setTxtProgressBar(
          pb,
          progress_step
        )
      }
    }


    # -------------------------------------------------------------------------
    # 8. Summarize paired simulations
    # -------------------------------------------------------------------------

    successful <-
      !is.na(lambda_difference) &
      !is.na(same_conclusion)

    n_successful <-
      sum(successful)


    if (n_successful == 0) {

      consistency <-
        NA_character_

      signal_bias <-
        NA_real_

      signal_sd <-
        NA_real_

    } else {

      consistency <- paste0(
        round(
          mean(
            same_conclusion[
              successful
            ]
          ) * 100,
          2
        ),
        " %"
      )


      signal_bias <- mean(
        lambda_difference[
          successful
        ]
      )


      signal_sd <- if (n_successful > 1) {

        stats::sd(
          lambda_difference[
            successful
          ]
        )

      } else {

        NA_real_
      }
    }


    sensitivity_results[[i]] <- data.frame(
      simulation_lambda = simulation_lambda,
      consistency = consistency,
      signal_bias = signal_bias,
      signal_sd = signal_sd,
      n_successful = n_successful
    )
  }


  # ---------------------------------------------------------------------------
  # 9. Close progress bar
  # ---------------------------------------------------------------------------

  if (!is.null(pb)) {
    close(pb)
    pb <- NULL
  }


  # ---------------------------------------------------------------------------
  # 10. Return results
  # ---------------------------------------------------------------------------

  sensitivity_results <- do.call(
    rbind,
    sensitivity_results
  )


  results <- cbind(
    observed_results,
    sensitivity_results
  )


  rownames(results) <- NULL


  attr(
    results,
    "n_simulations"
  ) <- n_simulations

  attr(
    results,
    "alpha_level"
  ) <- alpha_level


  results
}
