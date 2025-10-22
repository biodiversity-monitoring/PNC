#' Test Robustness of Phylogenetic Niche Conservatism Analysis Across Multiple Communities
#'
#' This function evaluates the robustness of phylogenetic signal estimates across
#' multiple communities by simulating trait data with the same phylogenetic signal
#' strength as observed, applying the original missing data pattern, and testing
#' how consistently the statistical significance is recovered across multiple
#' simulations for each community.
#'
#' @param com A community matrix with sites as rows and species as columns
#' @param trait_data A data frame or matrix containing trait data with species as rows
#' @param phylo_tree A phylogenetic tree object of class "phylo"
#' @param methods Character string specifying method to use. Options: "lambda" or "K".
#'   Default is "lambda"
#' @param pca_axes Character vector specifying which PCA axes to include (e.g., c("PC1", "PC2")).
#'   Default is c("PC1", "PC2")
#' @param sig_levels Numeric vector of significance levels for marking results
#' @param min_abundance Minimum abundance threshold for including species
#' @param n_simulations Integer. Number of simulations to run for robustness testing.
#'   Default is 100
#' @param alpha_level Numeric. Significance level for statistical testing.
#'   Default is 0.05
#' @param tolerance Numeric. Acceptable difference between target and estimated signal
#'   values during trait simulation. Default is 0.05
#' @param verbose Logical indicating whether to show progress and warnings
#'
#' @return A data frame containing the original phylogenetic signal results with
#'   additional columns:
#'   \itemize{
#'     \item robustness: Percentage of simulations that maintain the same statistical
#'           significance conclusion as the original analysis
#'     \item signal_sd: Standard deviation of phylogenetic signal values across
#'           successful simulations
#'   }
#'
#' @examples
#' \dontrun{
#' # Load example data
#' data("HimalayanBirds")
#' str(HimalayanBirds)
#' data("AVONET")
#' head(AVONET)
#'
#' # species level
#' sp <- colnames(HimalayanBirds$com)
#' sp
#' subtraits <- extract_traits(sp, AVONET, rank = "species")
#' head(subtraits)
#' coverage(subtraits)
#' pnc(subtraits, HimalayanBirds$phy_species, methods = "lambda", pca_axes = c("PC1", "PC2"))
#' pnc_robustness(subtraits, HimalayanBirds$phy_species,
#'                methods = "lambda", pca_axes = NULL, n_simulations = 100)
#'
#' compnc(com = HimalayanBirds$com, subtraits, HimalayanBirds$phy_species,
#'        methods = "lambda", pca_axes = NULL)
#'
#' # Test robustness of phylogenetic signal analysis
#' robust_results <- compnc_robustness(HimalayanBirds$com,
#'                                     subtraits,
#'                                     HimalayanBirds$phy_species,
#'                                     methods = "lambda",
#'                                     pca_axes = NULL,
#'                                     n_simulations = 100)
#' robust_results
#' }
#'
#' @export
#' @importFrom phytools phylosig
#' @importFrom ape drop.tip
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom stats complete.cases sd prcomp
compnc_robustness <- function(com, trait_data, phylo_tree,
                              methods = "lambda",
                              pca_axes = c("PC1", "PC2"),
                              sig_levels = c(0.001, 0.01, 0.05),
                              min_abundance = 0,
                              n_simulations = 100,
                              alpha_level = 0.05,
                              tolerance = 0.05,
                              verbose = TRUE) {
  # Step 1: Baseline Calculation using compnc
  if (verbose) {
    cat("Step 1: Performing baseline phylogenetic niche conservatism analysis...\n")
  }
  baseline_results <- compnc(com = com, trait_data = trait_data, phylo_tree = phylo_tree,
                             methods = methods, pca_axes = pca_axes, sig_levels = sig_levels,
                             min_abundance = min_abundance, nsim = 1000, verbose = verbose)
  # Get unique combinations of plot and trait for simulation
  unique_combinations <- unique(baseline_results[, c("plot", "trait")])
  if (verbose) {
    cat("Step 2: Starting robustness testing for", nrow(unique_combinations),
        "plot-trait combinations...\n")
  }
  # Initialize progress bar
  total_simulations <- nrow(unique_combinations) * n_simulations
  pb <- utils::txtProgressBar(min = 0, max = total_simulations, style = 3)
  progress_counter <- 0
  # Storage for simulation results
  simulation_results <- list()
  # Handle plot names
  plot_names <- base::rownames(com)
  if (base::is.null(plot_names)) {
    plot_names <- base::paste0("plot", base::sprintf("%03d", 1:base::nrow(com)))
    base::rownames(com) <- plot_names
  }
  # Main simulation loop for each plot-trait combination
  for (combo_i in 1:nrow(unique_combinations)) {
    plot_name <- unique_combinations$plot[combo_i]
    trait_name <- unique_combinations$trait[combo_i]
    # Get baseline results for this combination
    baseline_row <- baseline_results[baseline_results$plot == plot_name &
                                       baseline_results$trait == trait_name, ]
    if (nrow(baseline_row) == 0) {
      progress_counter <- progress_counter + n_simulations
      utils::setTxtProgressBar(pb, progress_counter)
      next
    }
    signal_obs <- baseline_row$signal[1]
    p_obs <- baseline_row$p[1]
    # Skip if baseline calculation failed
    if (is.na(signal_obs)) {
      simulation_results[[paste(plot_name, trait_name, sep = "_")]] <- list(
        signal_obs = NA,
        p_obs = NA,
        p_sim_values = rep(NA, n_simulations),
        signal_sim_values = rep(NA, n_simulations),
        n_successful = 0
      )
      progress_counter <- progress_counter + n_simulations
      utils::setTxtProgressBar(pb, progress_counter)
      next
    }
    # Get community data for this plot
    plot_abundances <- com[plot_name, ]
    present_species <- base::names(plot_abundances)[plot_abundances > min_abundance]
    # Match species with trait data
    species_indices <- base::match(present_species, base::rownames(trait_data))
    valid_indices <- !base::is.na(species_indices)
    if (base::sum(valid_indices) < 4) {
      simulation_results[[paste(plot_name, trait_name, sep = "_")]] <- list(
        signal_obs = signal_obs,
        p_obs = p_obs,
        p_sim_values = rep(NA, n_simulations),
        signal_sim_values = rep(NA, n_simulations),
        n_successful = 0
      )
      progress_counter <- progress_counter + n_simulations
      utils::setTxtProgressBar(pb, progress_counter)
      next
    }
    # Extract trait data for community species
    plot_trait_data <- trait_data[species_indices[valid_indices], , drop = FALSE]
    valid_species <- present_species[valid_indices]
    base::rownames(plot_trait_data) <- valid_species
    # Find common species between traits and phylogeny
    tree_species <- phylo_tree$tip.label
    common_species <- base::intersect(valid_species, tree_species)
    if (base::length(common_species) < 4) {
      simulation_results[[paste(plot_name, trait_name, sep = "_")]] <- list(
        signal_obs = signal_obs,
        p_obs = p_obs,
        p_sim_values = rep(NA, n_simulations),
        signal_sim_values = rep(NA, n_simulations),
        n_successful = 0
      )
      progress_counter <- progress_counter + n_simulations
      utils::setTxtProgressBar(pb, progress_counter)
      next
    }
    # Prune phylogenetic tree
    pruned_tree <- base::tryCatch({
      ape::drop.tip(phylo_tree, phylo_tree$tip.label[!phylo_tree$tip.label %in% common_species])
    }, error = function(e) {
      return(NULL)
    })
    if (base::is.null(pruned_tree)) {
      simulation_results[[paste(plot_name, trait_name, sep = "_")]] <- list(
        signal_obs = signal_obs,
        p_obs = p_obs,
        p_sim_values = rep(NA, n_simulations),
        signal_sim_values = rep(NA, n_simulations),
        n_successful = 0
      )
      progress_counter <- progress_counter + n_simulations
      utils::setTxtProgressBar(pb, progress_counter)
      next
    }
    # Prepare analysis data
    analysis_trait_data <- plot_trait_data[common_species, , drop = FALSE]
    # Handle PCA vs original traits for missing data pattern
    if (trait_name %in% pca_axes) {
      # For PCA axes, use complete cases missing pattern
      complete_cases <- stats::complete.cases(analysis_trait_data)
      na_pattern <- !complete_cases
      original_species <- base::rownames(analysis_trait_data)[complete_cases]
    } else {
      # For original traits
      trait_values <- analysis_trait_data[[trait_name]]
      na_pattern <- base::is.na(trait_values)
      original_species <- base::rownames(analysis_trait_data)[!na_pattern]
    }
    # Storage for this combination's simulation results
    p_sim_values <- numeric(n_simulations)
    signal_sim_values <- numeric(n_simulations)
    # Run simulations for this plot-trait combination
    for (sim_i in 1:n_simulations) {
      progress_counter <- progress_counter + 1
      utils::setTxtProgressBar(pb, progress_counter)
      base::tryCatch({
        # Simulate ideal data with target phylogenetic signal
        if (methods == "lambda") {
          ideal_trait <- simulate_lambda_trait(signal_obs, phylo_tree, tolerance = tolerance)
        } else if (methods == "K") {
          ideal_trait <- simulate_K_trait(signal_obs, phylo_tree, tolerance = tolerance)
        }
        if (base::is.null(ideal_trait)) {
          p_sim_values[sim_i] <- NA
          signal_sim_values[sim_i] <- NA
          next
        }
        # Apply exact missing pattern to simulated data
        simulated_trait <- ideal_trait[base::rownames(analysis_trait_data), , drop = FALSE]
        simulated_trait[na_pattern, ] <- NA
        # Prepare simulated trait vector for analysis
        trait_vector <- simulated_trait[, 1]
        base::names(trait_vector) <- base::rownames(simulated_trait)
        # Data cleaning and preparation
        valid_indices <- !base::is.na(trait_vector)
        clean_trait <- trait_vector[valid_indices]
        valid_species <- base::names(clean_trait)
        final_common_species <- base::intersect(valid_species, pruned_tree$tip.label)
        if (base::length(final_common_species) < 4) {
          p_sim_values[sim_i] <- NA
          signal_sim_values[sim_i] <- NA
          next
        }
        # Final tree pruning
        final_tree <- ape::drop.tip(pruned_tree,
                                    pruned_tree$tip.label[!pruned_tree$tip.label %in% final_common_species])
        # Prepare final trait vector
        final_trait <- clean_trait[final_common_species]
        final_trait <- final_trait[final_tree$tip.label]
        # Calculate phylogenetic signal
        if (methods == "lambda") {
          result <- phytools::phylosig(final_tree, final_trait,
                                       method = "lambda", test = TRUE, nsim = 1000)
          signal_sim_values[sim_i] <- result$lambda
        } else if (methods == "K") {
          result <- phytools::phylosig(final_tree, final_trait,
                                       method = "K", test = TRUE, nsim = 1000)
          signal_sim_values[sim_i] <- result$K
        }
        p_sim_values[sim_i] <- result$P
      }, error = function(e) {
        p_sim_values[sim_i] <- NA
        signal_sim_values[sim_i] <- NA
      })
    }
    # Store simulation results for this plot-trait combination
    simulation_results[[paste(plot_name, trait_name, sep = "_")]] <- list(
      signal_obs = signal_obs,
      p_obs = p_obs,
      p_sim_values = p_sim_values,
      signal_sim_values = signal_sim_values,
      n_successful = base::sum(!base::is.na(p_sim_values))
    )
  }
  base::close(pb)
  if (verbose) {
    cat("\nStep 3: Calculating robustness statistics...\n")
  }
  # Step 3: Calculate robustness statistics
  robustness_column <- numeric(nrow(baseline_results))
  signal_sd_column <- numeric(nrow(baseline_results))
  for (i in 1:nrow(baseline_results)) {
    plot_name <- baseline_results$plot[i]
    trait_name <- baseline_results$trait[i]
    key <- paste(plot_name, trait_name, sep = "_")
    if (!key %in% names(simulation_results)) {
      robustness_column[i] <- NA
      signal_sd_column[i] <- NA
      next
    }
    result <- simulation_results[[key]]
    if (is.na(result$signal_obs)) {
      robustness_column[i] <- NA
      signal_sd_column[i] <- NA
      next
    }
    # Calculate signal standard deviation
    valid_signals <- result$signal_sim_values[!base::is.na(result$signal_sim_values)]
    if (base::length(valid_signals) > 1) {
      signal_sd_column[i] <- base::round(stats::sd(valid_signals), 4)
    } else {
      signal_sd_column[i] <- NA
    }
    # Calculate robustness percentage
    valid_p_values <- result$p_sim_values[!base::is.na(result$p_sim_values)]
    if (base::length(valid_p_values) == 0) {
      robustness_column[i] <- NA
      next
    }
    # Determine consistency based on significance
    if (result$p_obs < alpha_level) {
      n_consistent <- base::sum(valid_p_values < alpha_level)
    } else {
      n_consistent <- base::sum(valid_p_values >= alpha_level)
    }
    robustness <- n_consistent / base::length(valid_p_values)
    robustness_column[i] <- base::paste0(base::round(robustness * 100, 2), " %")
  }
  # Add robustness and signal SD columns to baseline results
  final_results <- baseline_results
  final_results$robustness <- robustness_column
  final_results$signal_sd <- signal_sd_column
  # Set additional attributes
  base::attr(final_results, "n_simulations") <- n_simulations
  base::attr(final_results, "alpha_level") <- alpha_level
  base::attr(final_results, "tolerance") <- tolerance
  if (verbose) {
    cat("Robustness analysis complete!\n")
  }
  return(final_results)
}
