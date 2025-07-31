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
#' @param com A community matrix where rows represent sampling plots/sites and columns represent
#'   species. Values indicate species abundance or presence/absence.
#' @param trait_data A data frame containing trait measurements for species, where rows are species
#'   and columns are traits. Species names should match those in the community matrix.
#' @param phylo_tree A phylogenetic tree object (class "phylo") containing the evolutionary
#'   relationships among species. Tip labels should match species names in the community matrix.
#' @param methods Character vector specifying which phylogenetic signal methods to use. Available
#'   options include "Lambda" (Pagel's lambda), "K" (Blomberg's K), "K.star" (standardized K),
#'   "I" (Moran's I), and "Cmean" (mean phylogenetic distance). Default is all methods.
#' @param pca_axes Character vector specifying which principal component axes to include in the
#'   analysis (e.g., "PC1", "PC2"). These are calculated automatically from the trait data.
#' @param sig_levels Numeric vector of significance levels for statistical testing. Default is
#'   c(0.001, 0.01, 0.05) corresponding to ***, **, and * significance markers.
#' @param min_abundance Minimum abundance threshold for including species in the analysis. Species
#'   with abundance below this value are excluded. Default is 0.
#' @param reps Number of randomizations for statistical testing. Default is 999.
#' @param verbose Logical value indicating whether to display progress information during analysis.
#'   Default is TRUE.
#'
#' @return A data frame with the following columns:
#'   \describe{
#'     \item{plot}{Name or identifier of the community/plot}
#'     \item{trait}{Name of the trait being analyzed}
#'     \item{coverage}{Percentage of species in the community with available trait data}
#'     \item{n_sp}{Number of species included in the phylogenetic signal analysis}
#'     \item{signal}{Calculated phylogenetic signal value}
#'     \item{p}{P-value from randomization test}
#'     \item{significance}{Significance level marker (*** for p ≤ 0.001, ** for p ≤ 0.01, * for p ≤ 0.05, ns for non-significant)}
#'     \item{method}{Statistical method used for calculation}
#'     \item{n_sp_in_plot}{Total number of species present in the community}
#'   }
#'
#' @details
#' The function performs the following steps for each community:
#' \enumerate{
#'   \item Species Filtering: Identifies species present in each community based on the minimum abundance threshold.
#'   \item Data Matching: Matches species between community data, trait data, and phylogenetic tree,
#'         ensuring only species with complete information are included.
#'   \item Tree Pruning: Removes species from the phylogenetic tree that are not present in the
#'         current community.
#'   \item PCA Analysis: Performs principal component analysis on trait data to extract major axes
#'         of variation.
#'   \item Phylogenetic Signal Calculation: Computes phylogenetic signal for each trait using the
#'         specified methods.
#'   \item Statistical Testing: Evaluates the significance of phylogenetic signal through
#'         randomization tests.
#' }
#'
#' The function requires a minimum of 4 species with complete data for meaningful analysis.
#' Results are stored with analysis parameters as attributes.
#'
#' @examples
#' # Load example data
#' data(BCI)
#' data(TRY)
#'
#' # Extract trait data
#' sp <- colnames(BCI$com)
#' subtraits <- extract_traits(sp, TRY, rank = "species",
#'                            traits = c("LA", "LMA", "LeafN", "PlantHeight", "SeedMass", "SSD"))
#'
#' compnc(com = BCI$com, subtraits, BCI$phy_species, methods = "Lambda", pca_axes = NULL)
#'
#' @export
compnc <- function(com, trait_data, phylo_tree,
                   methods = c("Lambda", "K", "K.star", "I", "Cmean"),
                   pca_axes = c("PC1", "PC2"),
                   sig_levels = c(0.001, 0.01, 0.05),
                   min_abundance = 0,
                   reps = 999,
                   verbose = TRUE) {
  if (!is.data.frame(trait_data) && !is.matrix(trait_data)) {
    stop("The trait_data must be a data frame or a matrix!")
  }
  #if (any(is.na(trait_data)) && verbose) {
  #  warning(paste0("Trait dataset contains NA values."))
  #}
  required_packages <- c("phylobase", "phylosignal", "ape")
  for (pkg in required_packages) {
    if (!base::requireNamespace(pkg, quietly = TRUE)) {
      base::stop(base::paste("Package", pkg, "is required but not installed."))
    }
  }
  if (!base::requireNamespace("progress", quietly = TRUE)) {
    base::stop("Package 'progress' is required but not installed. Please install it using: install.packages('progress')")
  }
  available_methods <- c("Lambda", "K", "K.star", "I", "Cmean")
  methods <- base::match.arg(methods, available_methods, several.ok = TRUE)
  plot_names <- base::rownames(com)
  if (base::is.null(plot_names)) {
    plot_names <- base::paste0("plot", base::sprintf("%03d", 1:base::nrow(com)))
    base::rownames(com) <- plot_names
  }
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
  available_axes <- base::paste0("PC", 1:base::ncol(trait_data))
  pca_axes <- base::intersect(pca_axes, available_axes)
  get_significance <- function(p_value, levels = sig_levels) {
    if (base::is.na(p_value)) return("")
    if (p_value <= levels[1]) return("***")
    else if (p_value <= levels[2]) return("**")
    else if (p_value <= levels[3]) return("*")
    else return("ns")
  }
  calculate_phylo_signal <- function(tree, trait_vector, method, reps) {
    base::tryCatch({
      p4d <- phylobase::phylo4d(tree, tip.data = base::data.frame(trait = trait_vector))
      method_map <- c("Lambda" = "Lambda", "K" = "K", "K.star" = "K.star",
                      "I" = "I", "Cmean" = "Cmean")
      phylo_method <- method_map[method]
      result <- phylosignal::phyloSignal(p4d, methods = phylo_method, reps = reps)
      signal_value <- result$stat[1, phylo_method]
      p_value <- result$pvalue[1, phylo_method]
      return(base::list(signal = signal_value, p = p_value))
    }, error = function(e) {
      if (verbose) {
        base::warning(base::paste("Error calculating", method, ":", e$message))
      }
      return(base::list(signal = NA, p = NA))
    })
  }
  total_plots <- base::nrow(com)
  successful_plots <- 0
  if (verbose) {
    base::cat("Analyzing", total_plots, "communities for phylogenetic niche conservatism...\n")
  }
  if (verbose) {
    pb <- progress::progress_bar$new(
      format = "  Analysis Progress [:bar] :percent :current/:total completed :elapsed remaining :eta",
      total = total_plots,
      clear = FALSE,
      width = 80
    )
  }
  for (i in 1:base::nrow(com)) {
    plot_name <- plot_names[i]
    plot_abundances <- com[i, ]
    present_species <- base::names(plot_abundances)[plot_abundances > min_abundance]

    if (base::length(present_species) < 4) {
      if (verbose) {
        pb$tick(tokens = base::list(current = i, total = total_plots))
      }
      next
    }
    species_indices <- base::match(present_species, base::rownames(trait_data))
    valid_indices <- !base::is.na(species_indices)
    if (base::sum(valid_indices) < 4) {
      if (verbose) {
        pb$tick(tokens = base::list(current = i, total = total_plots))
      }
      next
    }
    plot_trait_data <- trait_data[species_indices[valid_indices], , drop = FALSE]
    valid_species <- present_species[valid_indices]
    base::rownames(plot_trait_data) <- valid_species
    tree_species <- phylo_tree$tip.label
    common_species <- base::intersect(valid_species, tree_species)

    if (base::length(common_species) < 4) {
      if (verbose) {
        pb$tick(tokens = base::list(current = i, total = total_plots))
      }
      next
    }
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
        pb$tick(tokens = base::list(current = i, total = total_plots))
      }
      next
    }
    analysis_trait_data <- plot_trait_data[common_species, , drop = FALSE]
    trait_names <- base::colnames(analysis_trait_data)
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
    for (trait_name in trait_names) {
      if (trait_name %in% pca_axes) {
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
        trait_values <- analysis_trait_data[[trait_name]]
        valid_indices <- !base::is.na(trait_values)
        trait_values <- trait_values[valid_indices]
        valid_trait_species <- base::rownames(analysis_trait_data)[valid_indices]
        coverage <- base::paste0(base::round((base::length(trait_values) / base::nrow(analysis_trait_data)) * 100, 2), " %")
      }
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
      final_tree <- ape::drop.tip(pruned_tree,
                                  pruned_tree$tip.label[!pruned_tree$tip.label %in% final_common_species])
      if (trait_name %in% pca_axes) {
        final_trait <- trait_values[valid_trait_species %in% final_common_species]
        base::names(final_trait) <- valid_trait_species[valid_trait_species %in% final_common_species]
      } else {
        final_trait <- trait_values[valid_trait_species %in% final_common_species]
        base::names(final_trait) <- valid_trait_species[valid_trait_species %in% final_common_species]
      }
      final_trait <- final_trait[final_tree$tip.label]
      for (method in methods) {
        result <- calculate_phylo_signal(final_tree, final_trait, method, reps)
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
    all_results <- base::rbind(all_results, plot_results)
    successful_plots <- successful_plots + 1
    if (verbose) {
      pb$tick(tokens = base::list(current = i, total = total_plots))
    }
  }
  base::attr(all_results, "methods") <- methods
  base::attr(all_results, "pca_axes") <- pca_axes
  base::attr(all_results, "sig_levels") <- sig_levels
  base::attr(all_results, "min_abundance") <- min_abundance
  base::attr(all_results, "reps") <- reps
  base::attr(all_results, "total_plots") <- total_plots
  base::attr(all_results, "analyzed_plots") <- successful_plots
  if (verbose) {
    base::cat("\nAnalysis complete! Successfully analyzed", successful_plots, "communities out of", total_plots, "communities\n")
  }
  return(all_results)
}

