#' Assess Sensitivity of Phylogenetic Signal to Missing Trait Data
#' Across Multiple Communities
#'
#' `compnc_robustness()` evaluates how observed patterns of missing trait data
#' affect inference of phylogenetic signal across multiple communities based
#' on Pagel's lambda.
#'
#' For each community and trait, the observed lambda is first estimated using
#' `compnc()`. Complete trait data are then simulated under Brownian motion on
#' the local reference phylogeny transformed according to the observed lambda.
#' Each simulated trait is analyzed twice: first using the complete local
#' species pool represented in the phylogeny, and then after applying the
#' observed missing-data pattern of the focal trait.
#'
#' If the observed lambda exceeds the maximum value supported by the local
#' reference phylogeny, the largest valid lambda for that phylogeny is used
#' for simulation. The observed lambda estimate itself is not modified.
#'
#' The paired difference between the two analyses quantifies the effect of
#' trait-data incompleteness on phylogenetic signal estimation.
#'
#' @param com A community matrix or data frame with communities as rows and
#'   species as columns. Cell values represent species abundance or occurrence.
#' @param trait_data A data frame or matrix containing continuous trait data,
#'   with species as rows and traits as columns. Species names must be stored
#'   as row names.
#' @param phylo_tree A phylogenetic tree object of class `"phylo"`.
#' @param min_abundance Minimum abundance required for a species to be
#'   considered present in a community. Default is 0.
#' @param n_simulations Number of paired simulations performed for each
#'   community-trait combination. Default is 100.
#' @param alpha_level Significance threshold used to compare statistical
#'   conclusions between complete and incomplete simulated data.
#'   Default is 0.05.
#' @param verbose Logical. If `TRUE`, a progress bar and warnings are shown.
#'   Default is `TRUE`.
#'
#' @return A data frame containing the observed phylogenetic signal results
#'   from `compnc()` and the following sensitivity metrics:
#'   \itemize{
#'     \item `simulation_lambda`: lambda value used to generate simulated trait
#'       data. This equals the observed lambda unless the latter exceeds the
#'       maximum supported by the local reference phylogeny.
#'     \item `consistency`: percentage of successful simulations in which
#'       complete and incomplete data produced the same significance
#'       conclusion.
#'     \item `signal_bias`: mean paired difference between lambda estimated
#'       from incomplete and complete simulated data
#'       (`lambda_missing - lambda_complete`).
#'     \item `signal_sd`: standard deviation of the paired lambda differences.
#'     \item `n_successful`: number of simulations in which both complete and
#'       incomplete analyses were successfully estimated.
#'   }
#'
#' @details
#' The analysis is performed separately for each community and original trait.
#' For each community, the reference species pool consists of species present
#' in the community and represented in the phylogeny. Species without trait
#' values are retained in the complete simulated data and removed only when
#' the observed missing-data pattern is applied.
#'
#' The simulation is conditional on the lambda estimated from the observed
#' data. The local reference phylogeny is transformed according to this lambda,
#' and trait values are simulated under Brownian motion on the transformed
#' tree.
#'
#' The maximum valid lambda depends on the branch-length structure of the
#' local reference phylogeny. Consequently, an observed lambda estimated from
#' a trait-specific pruned tree may occasionally exceed the maximum value that
#' can be applied to the larger local reference phylogeny. In this case,
#' `compnc_robustness()` identifies the largest valid lambda numerically and
#' uses that value only for simulation. A warning is issued when this
#' adjustment occurs.
#'
#' Each simulated trait is analyzed before and after applying the observed
#' missing-data pattern. Because the same simulated trait is used in both
#' analyses, differences between them reflect the effect of trait-data
#' incompleteness.
#'
#' Simulated lambda estimates are not forced to equal the observed lambda.
#' Their natural stochastic variation is retained.
#'
#' Repeated Pagel's lambda estimation is accelerated by caching the
#' phylogeny-specific covariance structure and using an algebraically
#' equivalent spectral representation of the likelihood. This changes the
#' numerical implementation but not the Pagel's lambda likelihood,
#' likelihood-ratio test, or P value.
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
#' set.seed(123)
#'
#' compnc_robustness(
#'   HimalayanBirds$com,
#'   subtraits,
#'   HimalayanBirds$phy_species,
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
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @export
compnc_robustness <- function(com,
                              trait_data,
                              phylo_tree,
                              min_abundance = 0,
                              n_simulations = 100,
                              alpha_level = 0.05,
                              verbose = TRUE) {

  # ---------------------------------------------------------------------------
  # 1. Check inputs
  # ---------------------------------------------------------------------------

  if (!inherits(phylo_tree, "phylo")) {
    stop("`phylo_tree` must be an object of class 'phylo'.")
  }

  com <- as.data.frame(
    com,
    check.names = FALSE
  )

  trait_data <- as.data.frame(
    trait_data,
    check.names = FALSE
  )

  if (is.null(rownames(trait_data))) {
    stop("`trait_data` must have species names as row names.")
  }

  if (anyDuplicated(rownames(trait_data))) {
    stop("Species names in `trait_data` must be unique.")
  }

  if (
    length(min_abundance) != 1 ||
    !is.finite(min_abundance)
  ) {
    stop("`min_abundance` must be a finite numeric value.")
  }

  if (
    length(n_simulations) != 1 ||
    !is.finite(n_simulations) ||
    n_simulations < 1 ||
    n_simulations %% 1 != 0
  ) {
    stop("`n_simulations` must be a positive integer.")
  }

  if (
    length(alpha_level) != 1 ||
    !is.finite(alpha_level) ||
    alpha_level <= 0 ||
    alpha_level >= 1
  ) {
    stop("`alpha_level` must be between 0 and 1.")
  }

  if (
    length(verbose) != 1 ||
    !is.logical(verbose) ||
    is.na(verbose)
  ) {
    stop("`verbose` must be TRUE or FALSE.")
  }

  n_simulations <- as.integer(n_simulations)


  # ---------------------------------------------------------------------------
  # 2. Community names
  # ---------------------------------------------------------------------------

  plot_names <- rownames(com)

  if (is.null(plot_names)) {

    plot_names <- paste0(
      "plot",
      sprintf(
        "%03d",
        seq_len(nrow(com))
      )
    )

    rownames(com) <- plot_names
  }


  # ---------------------------------------------------------------------------
  # 3. Calculate observed phylogenetic signal
  # ---------------------------------------------------------------------------

  observed_results <- compnc(
    com = com,
    trait_data = trait_data,
    phylo_tree = phylo_tree,
    methods = "lambda",
    pca_axes = NULL,
    min_abundance = min_abundance,
    verbose = FALSE
  )


  # ---------------------------------------------------------------------------
  # 4. Helper functions
  # ---------------------------------------------------------------------------

  # Standard Pagel's lambda estimator used only as a fallback.
  estimate_lambda_standard <- function(tree, trait) {

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

    if (is.null(result)) {
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

    C <- ape::vcv.phylo(tree)

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
          stop("Trait values are invalid.")
        }

        n <- length(trait)

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
              log(lambda_eigenvalues)
            )

          log_likelihood <-
            -quadratic_form / (2 * sigma2) -
            n * log(2 * pi) / 2 -
            (
              n * log(sigma2) +
                logdet_C_lambda
            ) / 2

          as.numeric(log_likelihood)
        }


        # Same ten-interval optimization strategy as phylosig().
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
          function(x) x$objective,
          numeric(1)
        )

        best_fit <- fits[[
          which.max(likelihoods)
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
        prepared,
        trait
      )

      if (
        is.finite(fast_result$lambda) &&
        is.finite(fast_result$p)
      ) {
        return(fast_result)
      }
    }

    estimate_lambda_standard(
      tree,
      trait
    )
  }


  # Apply Pagel's lambda transformation.
  #
  # Warnings caused by deliberately probing an invalid lambda during the
  # upper-bound search are suppressed. Other warnings remain visible.
  rescale_lambda <- function(tree, lambda) {

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
      is.null(transformed_tree$edge.length) ||
      any(!is.finite(transformed_tree$edge.length)) ||
      any(transformed_tree$edge.length < 0)
    ) {
      return(NULL)
    }

    transformed_tree
  }


  # Find the largest lambda supported by a local reference phylogeny.
  find_valid_lambda <- function(tree,
                                lambda_upper,
                                n_iterations = 60L) {

    lambda_lower <- 1

    tree_lower <- rescale_lambda(
      tree,
      lambda_lower
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
        tree,
        lambda_mid
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


  # ---------------------------------------------------------------------------
  # 5. Cache local reference phylogenies
  # ---------------------------------------------------------------------------

  # The complete local reference species pool is identical across traits
  # within a community. Cache the corresponding tree and spectral
  # decomposition so they are calculated only once per community.

  reference_cache <- vector(
    "list",
    length(plot_names)
  )

  names(reference_cache) <-
    plot_names


  get_reference_data <- function(plot_name) {

    if (!is.null(reference_cache[[plot_name]])) {
      return(
        reference_cache[[plot_name]]
      )
    }


    plot_abundance <- unlist(
      com[
        plot_name,
        ,
        drop = TRUE
      ],
      use.names = TRUE
    )

    present <- !is.na(plot_abundance) &
      plot_abundance > min_abundance

    present_species <- names(
      plot_abundance
    )[present]

    reference_species <- intersect(
      present_species,
      phylo_tree$tip.label
    )


    if (length(reference_species) < 4) {

      reference_data <- list(
        valid = FALSE,
        species = reference_species,
        tree = NULL,
        prepared = NULL
      )

      reference_cache[[plot_name]] <<-
        reference_data

      return(
        reference_data
      )
    }


    reference_tree <- ape::drop.tip(
      phylo_tree,
      setdiff(
        phylo_tree$tip.label,
        reference_species
      )
    )

    reference_species <-
      reference_tree$tip.label


    prepared_reference <- tryCatch(

      prepare_lambda_fast(
        reference_tree
      ),

      error = function(e) {

        if (verbose) {
          warning(
            "Fast Pagel's lambda preparation failed for plot `",
            plot_name,
            "`. Standard `phytools::phylosig()` will be used ",
            "for the complete-data analysis.",
            call. = FALSE
          )
        }

        NULL
      }
    )


    reference_data <- list(
      valid = TRUE,
      species = reference_species,
      tree = reference_tree,
      prepared = prepared_reference
    )


    reference_cache[[plot_name]] <<-
      reference_data


    reference_data
  }


  # ---------------------------------------------------------------------------
  # 6. Set up progress bar
  # ---------------------------------------------------------------------------

  total_steps <-
    nrow(observed_results) *
    n_simulations

  progress_step <- 0L

  if (verbose) {

    pb <- utils::txtProgressBar(
      min = 0,
      max = total_steps,
      style = 3
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

    plot_name <-
      observed_results$plot[i]

    trait_name <-
      observed_results$trait[i]

    lambda_obs <-
      observed_results$signal[i]


    # -------------------------------------------------------------------------
    # Skip combinations for which observed lambda could not be estimated
    # -------------------------------------------------------------------------

    if (!is.finite(lambda_obs)) {

      sensitivity_results[[i]] <- data.frame(
        simulation_lambda = NA_real_,
        consistency = NA_character_,
        signal_bias = NA_real_,
        signal_sd = NA_real_,
        n_successful = 0L
      )

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
    # Retrieve complete local reference species pool
    # -------------------------------------------------------------------------

    reference_data <- get_reference_data(
      plot_name
    )


    if (!reference_data$valid) {

      sensitivity_results[[i]] <- data.frame(
        simulation_lambda = NA_real_,
        consistency = NA_character_,
        signal_bias = NA_real_,
        signal_sd = NA_real_,
        n_successful = 0L
      )

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


    reference_species <-
      reference_data$species

    reference_tree <-
      reference_data$tree

    prepared_reference <-
      reference_data$prepared


    # -------------------------------------------------------------------------
    # Observed missing-data pattern
    # -------------------------------------------------------------------------

    trait_rows <- match(
      reference_species,
      rownames(trait_data)
    )

    has_trait <- rep(
      FALSE,
      length(reference_species)
    )

    matched <-
      !is.na(trait_rows)

    has_trait[matched] <- !is.na(
      trait_data[
        trait_rows[matched],
        trait_name
      ]
    )

    observed_species <- reference_species[
      has_trait
    ]


    if (length(observed_species) < 4) {

      sensitivity_results[[i]] <- data.frame(
        simulation_lambda = NA_real_,
        consistency = NA_character_,
        signal_bias = NA_real_,
        signal_sd = NA_real_,
        n_successful = 0L
      )

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
    # Prepare fast engine for trait-specific observed tree
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
            "` in plot `",
            plot_name,
            "`. Standard `phytools::phylosig()` will be used ",
            "for the incomplete-data analysis.",
            call. = FALSE
          )
        }

        NULL
      }
    )


    # -------------------------------------------------------------------------
    # Construct lambda-transformed local reference phylogeny
    # -------------------------------------------------------------------------

    simulation_lambda <-
      lambda_obs

    simulation_tree <- rescale_lambda(
      reference_tree,
      simulation_lambda
    )


    # If the observed lambda is not valid for the larger local reference tree,
    # identify the largest lambda that can be applied to that tree.
    if (
      is.null(simulation_tree) &&
      lambda_obs > 1
    ) {

      lambda_limit <- find_valid_lambda(
        reference_tree,
        lambda_obs
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
          "` in plot `",
          plot_name,
          "` (",
          formatC(
            lambda_obs,
            format = "f",
            digits = 6
          ),
          ") exceeds the maximum supported by the local reference ",
          "phylogeny. Lambda = ",
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
    # If transformation still fails, sensitivity cannot be evaluated
    # -------------------------------------------------------------------------

    if (is.null(simulation_tree)) {

      sensitivity_results[[i]] <- data.frame(
        simulation_lambda = NA_real_,
        consistency = NA_character_,
        signal_bias = NA_real_,
        signal_sd = NA_real_,
        n_successful = 0L
      )

      progress_step <-
        progress_step +
        n_simulations


      if (verbose) {

        warning(
          "Could not construct a valid lambda-transformed local ",
          "reference phylogeny for trait `",
          trait_name,
          "` in plot `",
          plot_name,
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


        # ---------------------------------------------------------------------
        # Complete-data analysis
        # ---------------------------------------------------------------------

        complete_result <- estimate_lambda(
          tree = reference_tree,
          prepared = prepared_reference,
          trait = simulated_trait
        )


        # ---------------------------------------------------------------------
        # Apply observed missing-data pattern
        # ---------------------------------------------------------------------

        missing_trait <- simulated_trait[
          observed_tree$tip.label
        ]


        # ---------------------------------------------------------------------
        # Incomplete-data analysis
        # ---------------------------------------------------------------------

        missing_result <- estimate_lambda(
          tree = observed_tree,
          prepared = prepared_observed,
          trait = missing_trait
        )


        # ---------------------------------------------------------------------
        # Keep only successful paired analyses
        # ---------------------------------------------------------------------

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


      # -----------------------------------------------------------------------
      # Update progress bar
      # -----------------------------------------------------------------------

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


      if (n_successful > 1) {

        signal_sd <- stats::sd(
          lambda_difference[
            successful
          ]
        )

      } else {

        signal_sd <-
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

  if (verbose) {
    close(pb)
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
    "min_abundance"
  ) <- min_abundance


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
