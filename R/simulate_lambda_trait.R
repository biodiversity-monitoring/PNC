#' Simulate Trait Data with Target Phylogenetic Signal (Lambda)
#'
#' This function generates trait data that matches a specified phylogenetic signal
#' strength (Pagel's lambda) through iterative simulation and testing.
#'
#' @param target_lambda Numeric. The desired phylogenetic signal strength (lambda value).
#'   Should be between 0 and 1.
#'   - 0: No phylogenetic signal (star phylogeny)
#'   - 1: Full phylogenetic signal (Brownian motion)
#' @param tree An object of class "phylo". The phylogenetic tree for trait simulation.
#' @param max_attempts Integer. Maximum number of simulation attempts before giving up.
#'   Default is 100000.
#' @param tolerance Numeric. Acceptable difference between target and estimated lambda.
#'   Default is 0.02.
#'
#' @return A data.frame with one column named 'trait' containing the simulated trait values.
#'   Row names correspond to tip labels from the phylogenetic tree. Returns NULL if
#'   the target lambda cannot be achieved within the specified tolerance and attempts.
#'
#' @details
#' The function works by:
#'
#' 1. Transforming the phylogenetic tree according to the target lambda value using phytools::rescale()
#'
#' 2. Simulating trait data using phytools::fastBM() on the transformed tree
#'
#' 3. Estimating the phylogenetic signal using phytools::phylosig()
#'
#' 4. Repeating until the estimated lambda is within tolerance of the target
#'
#' Special cases:
#' - When target_lambda = 0: Sets internal branch lengths to 0, keeping only terminal branches
#' - When target_lambda = 1: Uses the original tree without transformation
#'
#' @examples
#' \dontrun{
#' # Generate a random tree
#' tree <- ape::rtree(50)
#'
#' # Simulate trait with strong phylogenetic signal
#' trait_data <- simulate_lambda_trait(0.8, tree)
#'
#' # Verify the phylogenetic signal
#' trait_vector <- setNames(trait_data$trait, rownames(trait_data))
#' phytools::phylosig(tree, trait_vector, method = "lambda", test = TRUE)
#' }
#'
#' @note
#' The function may take considerable time to converge for certain lambda values,
#' especially those close to intermediate values. Consider adjusting
#' the tolerance parameter if convergence is slow.
#'
#'
#' @importFrom ape Ntip
#' @importFrom phytools fastBM phylosig rescale
#'
#' @export
simulate_lambda_trait <- function(target_lambda, tree, max_attempts = 100000, tolerance = 0.02) {
  attempt <- 1
  while (attempt <= max_attempts) {
    # Transform tree based on target lambda
    if (target_lambda == 0) {
      transformed_tree <- tree
      transformed_tree$edge.length <- ifelse(
        transformed_tree$edge[, 2] <= ape::Ntip(tree),
        transformed_tree$edge.length,
        0
      )
    } else if (target_lambda == 1) {
      transformed_tree <- tree
    } else {
      transformed_tree <- phytools::rescale(tree, "lambda", target_lambda)
    }
    # Simulate trait data
    trait <- phytools::fastBM(transformed_tree)
    # Test phylogenetic signal
    result <- phytools::phylosig(tree, trait, method = "lambda", test = TRUE)
    estimated_lambda <- result$lambda
    # Check if within tolerance
    if (abs(estimated_lambda - target_lambda) <= tolerance) {
      trait_df <- data.frame(
        trait = as.numeric(trait),
        row.names = names(trait)
      )
      return(trait_df)
    }
    attempt <- attempt + 1
  }
  # If we reach here, we failed to achieve target lambda
  warning("Failed to achieve target lambda within tolerance after ",
          max_attempts, " attempts")
  return(NULL)
}
