#' Assess Sensitivity of Phylogenetic Signal to Missing Trait Data
#' Across Multiple Communities
#'
#' `compnc_robustness()` evaluates how observed patterns of missing trait data
#' affect inference of phylogenetic signal across multiple communities based
#' on Pagel's lambda.
#'
#' For each community and trait, the observed lambda is first estimated using
#' `compnc()`. Complete trait data are then simulated under Brownian motion on
#' a local reference phylogeny transformed according to the observed lambda.
#' Each simulated trait is analyzed twice: first using the complete local
#' reference species pool, and then after applying the observed missing-data
#' pattern of the focal trait.
#'
#' The paired comparison isolates the effect of trait-data incompleteness on
#' phylogenetic signal estimation.
#'
#' @param com A community matrix or data frame with communities as rows and
#'   species as columns. Cell values represent species abundance or occurrence.
#' @param trait_data A data frame or matrix containing continuous trait data,
#'   with species as rows and traits as columns. Species names must be stored
#'   as row names.
#' @param phylo_tree A phylogenetic tree object of class `"phylo"` with branch
#'   lengths.
#' @param min_abundance Minimum abundance required for a species to be
#'   considered present in a community. A species is considered present when
#'   its abundance is greater than `min_abundance`. Default is 0.
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
#' The analysis is performed separately for each community and original trait.
#' For a given community, the local reference species pool consists of species
#' that are present in the community and represented in both `trait_data` and
#' `phylo_tree`. Species without observed values for the focal trait are retained
#' in the complete simulated data and removed only when the empirical
#' missing-data pattern is applied.
#'
#' For each community-trait combination, simulation is conditional on the
#' Pagel's lambda estimated from the observed data. The local reference
#' phylogeny is transformed according to this lambda, and trait values are
#' simulated under Brownian motion on the transformed tree.
#'
#' The maximum valid lambda depends on the branch-length structure of the
#' local reference phylogeny. An observed lambda estimated from a
#' trait-specific pruned tree may occasionally exceed the maximum value that
#' can be applied to the larger local reference phylogeny. In this case, the
#' largest valid lambda for that phylogeny is used only for simulation. The
#' observed lambda estimate itself is not modified.
#'
#' Each simulated trait is analyzed before and after applying the observed
#' missing-data pattern. Because both analyses use the same simulated trait,
#' their paired difference reflects the effect of trait-data incompleteness
#' rather than stochastic differences among simulations.
#'
#' Simulated lambda estimates are not forced to equal the observed lambda.
#' Their natural stochastic variation is retained.
#'
#' Repeated Pagel's lambda estimation is accelerated by caching
#' phylogeny-specific covariance structures and using an algebraically
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
#' @export
#' @import geiger
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

  if (!is.data.frame(com) && !is.matrix(com)) {
    stop("`com` must be a data frame or matrix.")
  }

  if (!is.data.frame(trait_data) && !is.matrix(trait_data)) {
    stop("`trait_data` must be a data frame or matrix.")
  }

  if (!inherits(phylo_tree, "phylo")) {
    stop("`phylo_tree` must be an object of class `phylo`.")
  }

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
    length(min_abundance) != 1 ||
    !is.finite(min_abundance) ||
    min_abundance < 0
  ) {
    stop(
      "`min_abundance` must be a single non-negative finite number."
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
  # 2. Calculate observed phylogenetic signal
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
  # 3. Helper functions
  # ---------------------------------------------------------------------------

  # Standard Pagel's lambda estimator used as a fallback
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


  # Prepare phylogeny-specific quantities for repeated lambda estimation
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


  # Fast Pagel's lambda estimator using cached spectral decomposition
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

        # Match the interval-search strategy used by phylosig()
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


  # Use the fast estimator and fall back to phylosig() if necessary
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


  # Apply Pagel's lambda transformation
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

  # Find the largest lambda supported by a local reference phylogeny
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


  # Standard empty sensitivity result
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
  # 4. Cache local reference phylogenies
  # ---------------------------------------------------------------------------

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

    # Complete local reference pool:
    # present in the community, represented in the trait table,
    # and represented in the phylogeny.
    reference_species <- Reduce(
      intersect,
      list(
        present_species,
        rownames(trait_data),
        phylo_tree$tip.label
      )
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
  # 5. Set up progress bar
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
  # 6. Run paired simulations
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
    # Retrieve complete local reference species pool
    # -------------------------------------------------------------------------

    reference_data <- get_reference_data(
      plot_name
    )

    if (!reference_data$valid) {

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


    reference_species <-
      reference_data$species

    reference_tree <-
      reference_data$tree

    prepared_reference <-
      reference_data$prepared


    # -------------------------------------------------------------------------
    # Identify the empirical missing-data pattern
    # -------------------------------------------------------------------------

    trait_values <- trait_data[
      reference_species,
      trait_name
    ]

    observed_species <- reference_species[
      !is.na(trait_values)
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
      tree = reference_tree,
      lambda = simulation_lambda
    )


    # The observed lambda may exceed the maximum supported by the larger
    # local reference tree.
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
    # 7. Summarize paired simulations
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
  # 8. Close progress bar
  # ---------------------------------------------------------------------------

  if (!is.null(pb)) {
    close(pb)
    pb <- NULL
  }


  # ---------------------------------------------------------------------------
  # 9. Return results
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
