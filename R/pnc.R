#' Calculate Phylogenetic Niche Conservatism for Single Community Analysis
#'
#' This function performs in-depth phylogenetic niche conservatism analysis for individual
#' communities by quantifying phylogenetic signal in trait data using multiple statistical methods.
#' The function integrates trait data preprocessing, phylogenetic tree manipulation,
#' optional principal component analysis, and robust statistical testing to provide
#' detailed insights into evolutionary constraints on niche occupation patterns.
#'
#' @param trait_data A data frame containing trait data with species as rows and traits as columns.
#'   Row names should correspond to species names that match the phylogenetic tree tip labels.
#' @param phylo_tree A phylogenetic tree object of class "phylo".
#'   Tip labels should match the species names in trait_data row names.
#' @param methods A character vector specifying which phylogenetic signal methods to use.
#'   Available options are: "Lambda", "K", "K.star", "I", "Cmean".
#'   Default is "Lambda".
#' @param pca_axes A character vector specifying which PCA axes to include in the analysis.
#'   Default is c("PC1", "PC2"). Set to NULL or empty vector to skip PCA analysis.
#' @param sig_levels A numeric vector of significance levels for marking results.
#'   Default is c(0.001, 0.01, 0.05) corresponding to ***, **, and * respectively.
#' @param reps An integer specifying the number of permutations for significance testing.
#'   Default is 999.
#' @param verbose A logical value indicating whether to display progress information and warnings.
#'   Default is TRUE.
#'
#' @return A data frame containing the phylogenetic signal results with the following columns:
#'   \describe{
#'     \item{trait}{Character. Name of the trait analyzed}
#'     \item{coverage}{Character. Percentage of species with valid data for the trait}
#'     \item{n_sp}{Integer. Number of species included in the analysis}
#'     \item{signal}{Numeric. Phylogenetic signal value}
#'     \item{p}{Numeric. P-value from permutation test}
#'     \item{significance}{Character. Significance level markers (***,**,*,ns)}
#'     \item{method}{Character. Method used for calculating phylogenetic signal}
#'   }
#'
#'   The returned object also contains the following attributes:
#'   \describe{
#'     \item{methods}{Character vector of methods used}
#'     \item{pca_axes}{Character vector of PCA axes requested}
#'     \item{pca_failed}{Logical indicating if PCA analysis failed}
#'     \item{sig_levels}{Numeric vector of significance levels used}
#'     \item{reps}{Integer number of permutations used}
#'     \item{pca_results}{Data frame of PCA scores (if PCA was successful)}
#'   }
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Validates input parameters and required packages
#'   \item Conducts PCA analysis on complete trait data (if requested)
#'   \item For each trait (original + PCA axes):
#'     \itemize{
#'       \item Identifies species with valid trait data
#'       \item Matches species between trait data and phylogenetic tree
#'       \item Prunes the phylogenetic tree to include only common species
#'       \item Calculates phylogenetic signal using specified methods
#'       \item Performs permutation tests for significance
#'     }
#'   \item Compiles results with coverage statistics and significance markers
#' }
#'
#' Phylogenetic Signal Methods:
#' \describe{
#'   \item{Lambda}{Pagel's lambda - measures the degree to which trait evolution follows Brownian motion}
#'   \item{K}{Blomberg's K - compares trait variation to that expected under Brownian motion}
#'   \item{K.star}{Modified K statistic accounting for tree structure}
#'   \item{I}{Moran's I - measures spatial autocorrelation in trait values}
#'   \item{Cmean}{Mean of squared changes - measures evolutionary rate}
#' }
#'
#' Data Requirements:
#' \itemize{
#'   \item Minimum 4 species with valid data for analysis
#'   \item Species names must match between trait data and phylogenetic tree
#'   \item Trait data should be numeric
#' }
#'
#' PCA Analysis:
#' If PCA axes are requested, the function performs PCA on complete cases only.
#' PCA is conducted with scaling and centering. If PCA fails or insufficient
#' complete data is available, the requested PCA axes will show NA values.
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
#' # Calculate phylogenetic signal using Lambda method
#' pnc(subtraits, BCI$phy_species, methods = "Lambda")
#'
#' # Calculate without PCA analysis
#' pnc(subtraits, BCI$phy_species, methods = "Lambda", pca_axes = NULL)
#' }
#'
#' @references
#' Münkemüller, T., Lavergne, S., Bzeznik, B., Dray, S., Jombart, T., Schiffers, K. and Thuiller, W. (2012). How to measure and test phylogenetic signal. Methods in Ecology and Evolution, 3(4), 743-756. \url{https://doi.org/10.1111/j.2041-210X.2012.00196.x}
#'
#' Keck, F., Rimet, F., Bouchez, A., & Franc, A. (2016). phylosignal: an R package to measure, test, and explore the phylogenetic signal. Ecology and Evolution, 6(9), 2774-2780. \url{https://doi.org/10.1002/ece3.2051}
#'
#' @export
#' @importFrom phylobase phylo4d
#' @importFrom phylosignal phyloSignal
#' @importFrom ape drop.tip
#' @importFrom stats prcomp complete.cases
#' @importFrom utils txtProgressBar setTxtProgressBar
pnc <- function(trait_data, phylo_tree,
                methods = "Lambda",
                pca_axes = c("PC1", "PC2"),
                sig_levels = c(0.001, 0.01, 0.05),
                reps = 999,
                verbose = TRUE) {
  required_packages <- c("phylobase", "phylosignal", "ape")
  for (pkg in required_packages) {
    if (!base::requireNamespace(pkg, quietly = TRUE)) {
      base::stop(base::paste("Package", pkg, "is required but not installed."))
    }
  }
  available_methods <- c("Lambda", "K", "K.star", "I", "Cmean")
  methods <- base::match.arg(methods, available_methods, several.ok = TRUE)
  available_axes <- base::paste0("PC", 1:base::ncol(trait_data))
  pca_axes <- base::intersect(pca_axes, available_axes)
  results <- base::data.frame(
    trait = base::character(),
    coverage = base::character(),
    n_sp = base::integer(),
    signal = base::numeric(),
    p = base::numeric(),
    significance = base::character(),
    method = base::character(),
    stringsAsFactors = FALSE
  )
  trait_names <- base::colnames(trait_data)
  n_traits <- base::length(trait_names)
  pca_results <- NULL
  pca_failed <- FALSE
  if (base::length(pca_axes) > 0) {
    complete_cases <- stats::complete.cases(trait_data)
    if (base::sum(complete_cases) >= 4) {
      complete_data <- trait_data[complete_cases, ]
      base::tryCatch({
        pca <- stats::prcomp(complete_data, scale. = TRUE, center = TRUE)
        pca_scores <- base::as.data.frame(pca$x)
        base::rownames(pca_scores) <- base::rownames(complete_data)
        selected_pca <- pca_scores[, pca_axes, drop = FALSE]
        pca_results <- selected_pca
        trait_names <- c(trait_names, base::colnames(selected_pca))
        n_traits <- base::length(trait_names)
      }, error = function(e) {
        if (verbose) {
          base::warning("PCA analysis failed: ", e$message)
        }
        pca_failed <<- TRUE
        trait_names <<- c(trait_names, pca_axes)
        n_traits <<- base::length(trait_names)
      })
    } else {
      if (verbose) {
        base::warning("Insufficient complete data samples for PCA analysis")
      }
      pca_failed <- TRUE
      trait_names <- c(trait_names, pca_axes)
      n_traits <- base::length(trait_names)
    }
  }
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
      method_map <- c("Lambda" = "Lambda", "K" = "K", "K.star" = "K.star", "I" = "I", "Cmean" = "Cmean")
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
  total_iterations <- n_traits * base::length(methods)
  if (verbose) {
    pb <- utils::txtProgressBar(min = 0, max = total_iterations, style = 3)
    current_iteration <- 0
  }
  for (i in base::seq_along(trait_names)) {
    trait_name <- trait_names[i]
    if (trait_name %in% pca_axes) {
      if (pca_failed || base::is.null(pca_results)) {
        for (method in methods) {
          results <- base::rbind(results, base::data.frame(
            trait = trait_name,
            coverage = NA,
            n_sp = NA,
            signal = NA,
            p = NA,
            significance = "",
            method = method
          ))
          if (verbose) {
            current_iteration <- current_iteration + 1
            utils::setTxtProgressBar(pb, current_iteration)
          }
        }
        next
      } else {
        trait_values <- pca_results[[trait_name]]
        valid_species <- base::rownames(pca_results)
        valid_trait <- trait_values
        coverage <- base::paste0(base::round((base::length(valid_trait) / base::nrow(trait_data)) * 100, 2), " %")
      }
    } else {
      trait_values <- trait_data[[trait_name]]
      valid_indices <- !base::is.na(trait_values)
      valid_trait <- trait_values[valid_indices]
      valid_species <- base::rownames(trait_data)[valid_indices]
      coverage <- base::paste0(base::round((base::length(valid_trait) / base::nrow(trait_data)) * 100, 2), " %")
    }
    n_sp <- base::length(valid_trait)
    if (n_sp < 4) {
      for (method in methods) {
        results <- base::rbind(results, base::data.frame(
          trait = trait_name,
          coverage = coverage,
          n_sp = n_sp,
          signal = NA,
          p = NA,
          significance = "",
          method = method
        ))
        if (verbose) {
          current_iteration <- current_iteration + 1
          utils::setTxtProgressBar(pb, current_iteration)
        }
      }
      next
    }
    tree_species <- phylo_tree$tip.label
    common_species <- base::intersect(valid_species, tree_species)
    if (base::length(common_species) < 4) {
      for (method in methods) {
        results <- base::rbind(results, base::data.frame(
          trait = trait_name,
          coverage = coverage,
          n_sp = n_sp,
          signal = NA,
          p = NA,
          significance = "",
          method = method
        ))
        if (verbose) {
          current_iteration <- current_iteration + 1
          utils::setTxtProgressBar(pb, current_iteration)
        }
      }
      next
    }
    pruned_tree <- ape::drop.tip(phylo_tree,
                                 phylo_tree$tip.label[!phylo_tree$tip.label %in% common_species])
    trait_for_analysis <- valid_trait[valid_species %in% common_species]
    base::names(trait_for_analysis) <- valid_species[valid_species %in% common_species]
    trait_for_analysis <- trait_for_analysis[pruned_tree$tip.label]
    for (method in methods) {
      result <- calculate_phylo_signal(pruned_tree, trait_for_analysis, method, reps)
      sig_mark <- get_significance(result$p)
      results <- base::rbind(results, base::data.frame(
        trait = trait_name,
        coverage = coverage,
        n_sp = base::length(common_species),
        signal = result$signal,
        p = result$p,
        significance = sig_mark,
        method = method
      ))
      if (verbose) {
        current_iteration <- current_iteration + 1
        utils::setTxtProgressBar(pb, current_iteration)
      }
    }
  }
  if (verbose) {
    base::close(pb)
  }
  base::attr(results, "methods") <- methods
  base::attr(results, "pca_axes") <- pca_axes
  base::attr(results, "pca_failed") <- pca_failed
  base::attr(results, "sig_levels") <- sig_levels
  base::attr(results, "reps") <- reps
  if (!base::is.null(pca_results)) {
    base::attr(results, "pca_results") <- pca_results
  }
  return(results)
}
