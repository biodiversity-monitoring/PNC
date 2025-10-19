#' Test Robustness of Phylogenetic Niche Conservatism Analysis
#'
#' This function evaluates the robustness of phylogenetic signal estimates by simulating
#' trait data with the same phylogenetic signal strength as observed, applying the
#' original missing data pattern, and testing how consistently the statistical
#' significance is recovered across multiple simulations.
#'
#' @param trait_data A data frame or matrix containing trait data with species as rows
#' @param phylo_tree A phylogenetic tree object of class "phylo"
#' @param methods Character string specifying method to use. Options: "lambda" or "K".
#'   Default is "lambda"
#' @param pca_axes Character vector specifying which PCA axes to include (e.g., c("PC1", "PC2")).
#'   Default is c("PC1", "PC2")
#' @param n_simulations Integer. Number of simulations to run for robustness testing.
#'   Default is 100
#' @param alpha_level Numeric. Significance level for statistical testing.
#'   Default is 0.05
#' @param tolerance Numeric. Acceptable difference between target and estimated signal
#'   values during trait simulation. Default is 0.05
#'
#' @return A data frame containing the original phylogenetic signal results with
#'   additional columns:
#'   \itemize{
#'     \item robustness: Percentage of simulations that maintain the same statistical
#'           significance conclusion as the original analysis
#'     \item signal_sd: Standard deviation of phylogenetic signal values across
#'           successful simulations
#'   }
#'   Returns the enhanced results from the baseline pnc() analysis
#'
#' @details
#' The robustness testing procedure involves:
#'
#' 1. Performing baseline phylogenetic signal analysis using pnc()
#'
#' 2. For each trait, simulating new trait data with the same phylogenetic signal
#'    strength as observed in the original data
#'
#' 3. Applying the exact missing data pattern from the original dataset to the
#'    simulated data
#'
#' 4. Re-testing phylogenetic signal on the simulated data and recording p-values
#'
#' 5. Calculating the percentage of simulations that maintain the same statistical
#'    significance conclusion (significant vs. non-significant)
#'
#' The function uses simulate_lambda_trait() or simulate_K_trait() internally to
#' generate trait data with target phylogenetic signal values.
#'
#' For PCA axes, the missing data pattern corresponds to complete cases from the
#' original trait matrix. For individual traits, the original missing pattern is
#' preserved exactly.
#'
#' @examples
#' \dontrun{
#' # Load example data
#' data(BCI)
#' data(TRY)
#'
#' # Extract trait data
#' sp <- colnames(BCI$com)
#' subtraits <- extract_traits(sp, TRY, rank = "species",
#'                           traits = c("LA", "LMA", "LeafN", "PlantHeight"))
#'
#' # Test robustness of phylogenetic signal analysis
#' robust_results <- pnc_robustness(subtraits, BCI$phy_species,
#'                                 methods = "lambda", n_simulations = 100)
#' robust_results
#' }
#'
#'
#' @importFrom phytools phylosig
#' @importFrom ape drop.tip
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom stats complete.cases sd
#'
#' @export
pnc_robustness <- function(trait_data,
                           phylo_tree,
                           methods = "lambda",
                           pca_axes = c("PC1", "PC2"),
                           n_simulations = 100,
                           alpha_level = 0.05,
                           tolerance = 0.05) {
  # Step 1: Baseline Calculation
  baseline_results <- pnc(trait_data, phylo_tree, methods = methods, pca_axes = pca_axes)
  # Get all trait names
  trait_names <- unique(baseline_results$trait)
  # Step 2: Simulation Analysis with progress bar
  pb <- txtProgressBar(min = 0, max = length(trait_names) * n_simulations, style = 3)
  progress_counter <- 0
  # Create storage for simulation results
  simulation_results <- list()
  # Start simulation for each trait
  for (trait_name in trait_names) {
    # Get baseline signal value for this trait
    baseline_row <- baseline_results[baseline_results$trait == trait_name, ]
    signal_obs <- baseline_row$signal
    p_obs <- baseline_row$p

    # Skip if baseline calculation failed
    if (is.na(signal_obs)) {
      simulation_results[[trait_name]] <- list(
        signal_obs = NA,
        p_obs = NA,
        p_sim_values = rep(NA, n_simulations),
        signal_sim_values = rep(NA, n_simulations),
        n_successful = 0
      )
      progress_counter <- progress_counter + n_simulations
      setTxtProgressBar(pb, progress_counter)
      next
    }
    # Get original data missing pattern
    if (trait_name %in% pca_axes) {
      # For PCA axes, use complete cases missing pattern
      complete_cases <- complete.cases(trait_data)
      original_species <- rownames(trait_data)[complete_cases]
      na_pattern <- !complete_cases
    } else {
      # For original traits
      trait_values <- trait_data[[trait_name]]
      na_pattern <- is.na(trait_values)
      original_species <- rownames(trait_data)[!na_pattern]
    }
    # Store p-values and signal values from each simulation
    p_sim_values <- numeric(n_simulations)
    signal_sim_values <- numeric(n_simulations)
    # Start simulation loop
    for (sim_i in 1:n_simulations) {
      progress_counter <- progress_counter + 1
      setTxtProgressBar(pb, progress_counter)
      tryCatch({
        # Step 2a: Simulate ideal data
        if (methods == "lambda") {
          ideal_trait <- simulate_lambda_trait(signal_obs, phylo_tree, tolerance = tolerance)
        } else if (methods == "K") {
          ideal_trait <- simulate_K_trait(signal_obs, phylo_tree, tolerance = tolerance)
        }
        if (is.null(ideal_trait)) {
          p_sim_values[sim_i] <- NA
          signal_sim_values[sim_i] <- NA
          next
        }
        # Step 2b: Apply exact missing pattern
        simulated_trait <- ideal_trait[rownames(trait_data), , drop = FALSE]
        simulated_trait[na_pattern, ] <- NA
        # Step 2c: Re-test phylogenetic signal of simulated data
        trait_vector <- simulated_trait[, 1]
        names(trait_vector) <- rownames(simulated_trait)
        # Data cleaning and preparation
        valid_indices <- !is.na(trait_vector)
        clean_trait <- trait_vector[valid_indices]
        valid_species <- names(clean_trait)
        tree_species <- phylo_tree$tip.label
        common_species <- intersect(valid_species, tree_species)
        if (length(common_species) < 4) {
          p_sim_values[sim_i] <- NA
          signal_sim_values[sim_i] <- NA
          next
        }
        # Prune phylogenetic tree
        pruned_tree <- ape::drop.tip(phylo_tree,
                                     phylo_tree$tip.label[!phylo_tree$tip.label %in% common_species])
        # Prepare trait data for analysis
        trait_for_analysis <- clean_trait[common_species]
        trait_for_analysis <- trait_for_analysis[pruned_tree$tip.label]
        # Calculate phylogenetic signal
        if (methods == "lambda") {
          result <- phytools::phylosig(pruned_tree, trait_for_analysis,
                                       method = "lambda", test = TRUE, nsim = 1000)
          signal_sim_values[sim_i] <- result$lambda
        } else if (methods == "K") {
          result <- phytools::phylosig(pruned_tree, trait_for_analysis,
                                       method = "K", test = TRUE, nsim = 1000)
          signal_sim_values[sim_i] <- result$K
        }
        p_sim_values[sim_i] <- result$P
      }, error = function(e) {
        p_sim_values[sim_i] <- NA
        signal_sim_values[sim_i] <- NA
      })
    }
    # Store simulation results for this trait
    simulation_results[[trait_name]] <- list(
      signal_obs = signal_obs,
      p_obs = p_obs,
      p_sim_values = p_sim_values,
      signal_sim_values = signal_sim_values,
      n_successful = sum(!is.na(p_sim_values))
    )
  }
  close(pb)
  # Step 3: Robustness Calculation
  # Add robustness and signal SD columns to baseline results
  robustness_column <- numeric(nrow(baseline_results))
  signal_sd_column <- numeric(nrow(baseline_results))
  for (i in 1:nrow(baseline_results)) {
    trait_name <- baseline_results$trait[i]
    result <- simulation_results[[trait_name]]
    if (is.na(result$signal_obs)) {
      robustness_column[i] <- NA
      signal_sd_column[i] <- NA
      next
    }
    valid_signals <- result$signal_sim_values[!is.na(result$signal_sim_values)]
    if (length(valid_signals) > 1) {
      signal_sd_column[i] <- round(sd(valid_signals), 4)
    } else {
      signal_sd_column[i] <- NA
    }
    valid_p_values <- result$p_sim_values[!is.na(result$p_sim_values)]
    if (length(valid_p_values) == 0) {
      robustness_column[i] <- NA
      next
    }
    if (result$p_obs < alpha_level) {
      n_consistent <- sum(valid_p_values < alpha_level)
    } else {
      n_consistent <- sum(valid_p_values >= alpha_level)
    }
    robustness <- n_consistent / length(valid_p_values)
    robustness_column[i] <- paste0(round(robustness * 100, 2), " %")
  }
  # Add robustness and signal SD columns to baseline results
  final_results <- baseline_results
  final_results$robustness <- robustness_column
  final_results$signal_sd <- signal_sd_column
  return(final_results)
}
