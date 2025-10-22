#' Analyze Phylogenetic Niche Conservatism in Ecological Communities
#'
#' This function performs in-depth phylogenetic niche conservatism analysis for
#' communities by quantifying phylogenetic signal in trait data using multiple statistical methods.
#' The function integrates trait data preprocessing, phylogenetic tree manipulation,
#' optional principal component analysis, and robust statistical testing to provide detailed insights
#' into evolutionary constraints on trait evolution.
#'
#' @param trait_data A data frame or matrix containing trait data with species as rows
#' @param phylo_tree A phylogenetic tree object of class "phylo"
#' @param methods Character vector specifying methods to use. Options: "lambda", "K"
#' @param pca_axes Character vector specifying which PCA axes to include (e.g., c("PC1", "PC2"))
#' @param sig_levels Numeric vector of significance levels for marking results
#' @param nsim Number of permutations for significance testing
#' @param verbose Logical indicating whether to show progress and warnings
#'
#' @return A data frame containing phylogenetic signal results
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
#'                             traits = c("LA", "LMA", "LeafN", "PlantHeight", "SeedMass", "SSD"))
#'
#' # Calculate phylogenetic signal using Lambda method
#' pnc(subtraits, BCI$phy_species, methods = "lambda")
#'
#' # Calculate without PCA analysis
#' pnc(subtraits, BCI$phy_species, methods = "lambda", pca_axes = NULL)
#' }
#'
#' @references
#' Münkemüller, T., Lavergne, S., Bzeznik, B., Dray, S., Jombart, T., Schiffers, K. and Thuiller, W. (2012).
#' How to measure and test phylogenetic signal. Methods in Ecology and Evolution, 3(4), 743-756.
#' \url{https://doi.org/10.1111/j.2041-210X.2012.00196.x}
#'
#' @export
#' @importFrom phytools phylosig
#' @importFrom ape drop.tip
#' @importFrom stats prcomp complete.cases
#' @importFrom utils txtProgressBar setTxtProgressBar
pnc <- function(trait_data, phylo_tree,
                methods = "lambda",
                pca_axes = c("PC1", "PC2"),
                sig_levels = c(0.001, 0.01, 0.05),
                nsim = 1000,
                verbose = TRUE) {
  if (!is.data.frame(trait_data) && !is.matrix(trait_data)) {
    stop("The trait_data must be a data frame or a matrix!")
  }
  required_packages <- c("phytools", "ape")
  for (pkg in required_packages) {
    if (!base::requireNamespace(pkg, quietly = TRUE)) {
      base::stop(base::paste("Package", pkg, "is required but not installed."))
    }
  }
  available_methods <- c("lambda", "K")
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
      result <- calculate_phylo_signal(pruned_tree, trait_for_analysis, method, nsim)
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
  base::attr(results, "nsim") <- nsim
  if (!base::is.null(pca_results)) {
    base::attr(results, "pca_results") <- pca_results
  }
  return(results)
}
