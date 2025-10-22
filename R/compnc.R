#' Calculate Phylogenetic Niche Conservatism Across Multiple Communities
#'
#' This function conducts comprehensive phylogenetic niche conservatism analysis across multiple
#' communities simultaneously. It evaluates phylogenetic signal for trait data across different
#' community assemblages using various statistical methods, enabling comparative assessment of
#' niche conservatism patterns among communities. The function processes community composition
#' matrices, species trait information, and phylogenetic trees to determine whether closely
#' related species consistently occupy similar ecological niches across different habitats
#' or sampling locations.
#'
#' @param com A community matrix with sites as rows and species as columns
#' @param trait_data A data frame or matrix containing trait data with species as rows
#' @param phylo_tree A phylogenetic tree object of class "phylo"
#' @param methods Character vector specifying methods to use. Options: "lambda", "K"
#' @param pca_axes Character vector specifying which PCA axes to include (e.g., c("PC1", "PC2"))
#' @param sig_levels Numeric vector of significance levels for marking results
#' @param min_abundance Minimum abundance threshold for including species
#' @param nsim Number of permutations for significance testing
#' @param verbose Logical indicating whether to show progress and warnings
#'
#' @return A data frame containing phylogenetic signal results for all communities
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
#'                            traits = c("LA", "LMA", "LeafN", "PlantHeight", "SeedMass", "SSD"))
#'
#' compnc(com = BCI$com, subtraits, BCI$phy_species, methods = "lambda", pca_axes = NULL)
#'}
#'
#' @export
#' @importFrom phytools phylosig
#' @importFrom ape drop.tip
#' @importFrom stats prcomp complete.cases
#' @importFrom utils txtProgressBar setTxtProgressBar
compnc <- function(com, trait_data, phylo_tree,
                   methods = c("lambda", "K"),
                   pca_axes = c("PC1", "PC2"),
                   sig_levels = c(0.001, 0.01, 0.05),
                   min_abundance = 0,
                   nsim = 1000,
                   verbose = TRUE) {
  # Input validation
  if (!is.data.frame(trait_data) && !is.matrix(trait_data)) {
    stop("The trait_data must be a data frame or a matrix!")
  }
  # Validate methods
  available_methods <- c("lambda", "K")
  methods <- base::match.arg(methods, available_methods, several.ok = TRUE)
  # Handle plot names
  plot_names <- base::rownames(com)
  if (base::is.null(plot_names)) {
    plot_names <- base::paste0("plot", base::sprintf("%03d", 1:base::nrow(com)))
    base::rownames(com) <- plot_names
  }
  # Initialize results data frame
  all_results <- base::data.frame(
    plot = base::character(),
    trait = base::character(),
    coverage = base::character(),
    n_sp = base::integer(),
    signal = base::numeric(),
    p = base::numeric(),
    significance = base::character(),
    method = base::character(),
    n_sp_in_plot = base::integer(),
    stringsAsFactors = FALSE
  )
  # Validate PCA axes
  available_axes <- base::paste0("PC", 1:base::ncol(trait_data))
  pca_axes <- base::intersect(pca_axes, available_axes)
  # Helper function for significance marking
  get_significance <- function(p_value, levels = sig_levels) {
    if (base::is.na(p_value)) return("")
    if (p_value <= levels[1]) return("***")
    else if (p_value <= levels[2]) return("**")
    else if (p_value <= levels[3]) return("*")
    else return("ns")
  }
  # Helper function for phylogenetic signal calculation using phytools
  calculate_phylo_signal <- function(tree, trait_vector, method, nsim) {
    base::tryCatch({
      if (method == "lambda") {
        result <- phytools::phylosig(tree, trait_vector, method = "lambda", test = TRUE, nsim = nsim)
        return(base::list(signal = result$lambda, p = result$P))
      } else if (method == "K") {
        result <- phytools::phylosig(tree, trait_vector, method = "K", test = TRUE, nsim = nsim)
        return(base::list(signal = result$K, p = result$P))
      }
    }, error = function(e) {
      if (verbose) {
        base::warning(base::paste("Error calculating", method, ":", e$message))
      }
      return(base::list(signal = NA, p = NA))
    })
  }
  # Progress tracking
  total_plots <- base::nrow(com)
  successful_plots <- 0
  if (verbose) {
    base::cat("Analyzing", total_plots, "communities for phylogenetic niche conservatism...\n")
    pb <- utils::txtProgressBar(min = 0, max = total_plots, style = 3)
  }
  # Main analysis loop for each community
  for (i in 1:base::nrow(com)) {
    plot_name <- plot_names[i]
    plot_abundances <- com[i, ]
    present_species <- base::names(plot_abundances)[plot_abundances > min_abundance]
    # Skip communities with too few species
    if (base::length(present_species) < 4) {
      if (verbose) {
        utils::setTxtProgressBar(pb, i)
      }
      next
    }
    # Match species with trait data
    species_indices <- base::match(present_species, base::rownames(trait_data))
    valid_indices <- !base::is.na(species_indices)
    if (base::sum(valid_indices) < 4) {
      if (verbose) {
        utils::setTxtProgressBar(pb, i)
      }
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
      if (verbose) {
        utils::setTxtProgressBar(pb, i)
      }
      next
    }
    # Prune phylogenetic tree
    pruned_tree <- base::tryCatch({
      ape::drop.tip(phylo_tree, phylo_tree$tip.label[!phylo_tree$tip.label %in% common_species])
    }, error = function(e) {
      if (verbose) {
        base::warning("Failed to prune tree for", plot_name, ":", e$message)
      }
      return(NULL)
    })
    if (base::is.null(pruned_tree)) {
      if (verbose) {
        utils::setTxtProgressBar(pb, i)
      }
      next
    }
    # Prepare analysis data
    analysis_trait_data <- plot_trait_data[common_species, , drop = FALSE]
    trait_names <- base::colnames(analysis_trait_data)
    # PCA analysis
    pca_results <- NULL
    pca_failed <- FALSE
    if (base::length(pca_axes) > 0) {
      complete_cases <- stats::complete.cases(analysis_trait_data)
      if (base::sum(complete_cases) >= 4) {
        complete_data <- analysis_trait_data[complete_cases, ]
        base::tryCatch({
          pca <- stats::prcomp(complete_data, scale. = TRUE, center = TRUE)
          pca_scores <- base::as.data.frame(pca$x)
          base::rownames(pca_scores) <- base::rownames(complete_data)
          selected_pca <- pca_scores[, pca_axes, drop = FALSE]
          pca_results <- selected_pca
          trait_names <- c(trait_names, base::colnames(selected_pca))
        }, error = function(e) {
          if (verbose) {
            base::warning("PCA analysis failed for", plot_name, ":", e$message)
          }
          pca_failed <- TRUE
        })
      } else {
        pca_failed <- TRUE
      }
    }
    # Initialize plot results
    plot_results <- base::data.frame(
      plot = base::character(),
      trait = base::character(),
      coverage = base::character(),
      n_sp = base::integer(),
      signal = base::numeric(),
      p = base::numeric(),
      significance = base::character(),
      method = base::character(),
      n_sp_in_plot = base::integer(),
      stringsAsFactors = FALSE
    )
    # Analyze each trait
    for (trait_name in trait_names) {
      if (trait_name %in% pca_axes) {
        # Handle PCA traits
        if (pca_failed || base::is.null(pca_results)) {
          for (method in methods) {
            plot_results <- base::rbind(plot_results, base::data.frame(
              plot = plot_name,
              trait = trait_name,
              coverage = NA,
              n_sp = NA,
              signal = NA,
              p = NA,
              significance = "",
              method = method,
              n_sp_in_plot = base::length(present_species)
            ))
          }
          next
        } else {
          trait_values <- pca_results[[trait_name]]
          valid_trait_species <- base::rownames(pca_results)
          coverage <- base::paste0(base::round((base::length(trait_values) / base::nrow(analysis_trait_data)) * 100, 2), " %")
        }
      } else {
        # Handle original traits
        trait_values <- analysis_trait_data[[trait_name]]
        valid_indices <- !base::is.na(trait_values)
        trait_values <- trait_values[valid_indices]
        valid_trait_species <- base::rownames(analysis_trait_data)[valid_indices]
        coverage <- base::paste0(base::round((base::length(trait_values) / base::nrow(analysis_trait_data)) * 100, 2), " %")
      }
      # Check minimum species requirement
      if (base::length(trait_values) < 4) {
        for (method in methods) {
          plot_results <- base::rbind(plot_results, base::data.frame(
            plot = plot_name,
            trait = trait_name,
            coverage = coverage,
            n_sp = base::length(trait_values),
            signal = NA,
            p = NA,
            significance = "",
            method = method,
            n_sp_in_plot = base::length(present_species)
          ))
        }
        next
      }
      # Final species matching with phylogeny
      final_common_species <- base::intersect(valid_trait_species, pruned_tree$tip.label)
      if (base::length(final_common_species) < 4) {
        for (method in methods) {
          plot_results <- base::rbind(plot_results, base::data.frame(
            plot = plot_name,
            trait = trait_name,
            coverage = coverage,
            n_sp = base::length(trait_values),
            signal = NA,
            p = NA,
            significance = "",
            method = method,
            n_sp_in_plot = base::length(present_species)
          ))
        }
        next
      }
      # Final tree pruning
      final_tree <- ape::drop.tip(pruned_tree,
                                  pruned_tree$tip.label[!pruned_tree$tip.label %in% final_common_species])

      # Prepare final trait vector
      final_trait <- trait_values[valid_trait_species %in% final_common_species]
      base::names(final_trait) <- valid_trait_species[valid_trait_species %in% final_common_species]
      final_trait <- final_trait[final_tree$tip.label]
      # Calculate phylogenetic signal for each method
      for (method in methods) {
        result <- calculate_phylo_signal(final_tree, final_trait, method, nsim)
        sig_mark <- get_significance(result$p)
        plot_results <- base::rbind(plot_results, base::data.frame(
          plot = plot_name,
          trait = trait_name,
          coverage = coverage,
          n_sp = base::length(final_common_species),
          signal = result$signal,
          p = result$p,
          significance = sig_mark,
          method = method,
          n_sp_in_plot = base::length(present_species)
        ))
      }
    }
    # Add plot results to overall results
    all_results <- base::rbind(all_results, plot_results)
    successful_plots <- successful_plots + 1
    if (verbose) {
      utils::setTxtProgressBar(pb, i)
    }
  }
  # Close progress bar
  if (verbose) {
    base::close(pb)
    base::cat("\nAnalysis complete! Successfully analyzed", successful_plots, "communities out of", total_plots, "communities\n")
  }
  # Set attributes
  base::attr(all_results, "methods") <- methods
  base::attr(all_results, "pca_axes") <- pca_axes
  base::attr(all_results, "sig_levels") <- sig_levels
  base::attr(all_results, "min_abundance") <- min_abundance
  base::attr(all_results, "nsim") <- nsim
  base::attr(all_results, "total_plots") <- total_plots
  base::attr(all_results, "analyzed_plots") <- successful_plots
  return(all_results)
}
