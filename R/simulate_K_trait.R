#' Simulate Trait Data with Target Phylogenetic Signal (Blomberg's K)
#'
#' This function generates trait data that matches a specified phylogenetic signal
#' strength (Blomberg's K) through iterative simulation and testing.
#'
#' @param target_K Numeric. The desired phylogenetic signal strength (K value).
#'   - K = 0: No phylogenetic signal (star phylogeny)
#'   - K = 1: Expected signal under Brownian motion evolution
#'   - K > 1: Stronger phylogenetic signal than expected under Brownian motion
#'   - 0 < K < 1: Weaker phylogenetic signal than expected under Brownian motion
#' @param tree An object of class "phylo". The phylogenetic tree for trait simulation.
#' @param max_attempts Integer. Maximum number of simulation attempts before giving up.
#'   Default is 100000.
#' @param tolerance Numeric. Acceptable difference between target and estimated K.
#'   Default is 0.02.
#'
#' @return A data.frame with one column named 'trait' containing the simulated trait values.
#'   Row names correspond to tip labels from the phylogenetic tree. Returns NULL if
#'   the target K cannot be achieved within the specified tolerance and attempts.
#'
#' @details
#' The function works by:
#'
#' 1. Transforming the phylogenetic tree according to the target K value
#'
#' 2. Simulating trait data using phytools::fastBM() on the transformed tree
#'
#' 3. Estimating the phylogenetic signal using phytools::phylosig()
#'
#' 4. Repeating until the estimated K is within tolerance of the target
#'
#' Tree transformation strategies:
#' - When target_K = 0: Creates a star phylogeny using ape::stree()
#' - When target_K = 1: Uses the original tree without transformation
#' - When target_K > 1: Scales all branch lengths by the target K value
#' - When 0 < target_K < 1: Interpolates between original tree and uniform branch lengths
#'
#' @examples
#' \dontrun{
#' # Generate a random tree
#' tree <- ape::rtree(50)
#'
#' # Simulate trait with expected Brownian motion signal
#' trait_data <- simulate_K_trait(0.9, tree)
#'
#' # Verify the phylogenetic signal
#' trait_vector <- setNames(trait_data$trait, rownames(trait_data))
#' phytools::phylosig(tree, trait_vector, method = "K", test = TRUE)
#' }
#'
#' @note
#' Blomberg's K measures the strength of phylogenetic signal relative to what would
#' be expected under a Brownian motion model of evolution. Unlike Pagel's lambda,
#' K can exceed 1, indicating stronger phylogenetic clustering than expected.
#'
#' The function may take considerable time to converge for certain K values.
#' Consider adjusting the tolerance parameter if convergence is slow.
#'
#' @importFrom ape Ntip stree
#' @importFrom phytools fastBM phylosig
#'
#' @export
simulate_K_trait <- function(target_K, tree, max_attempts = 100000, tolerance = 0.02) {
  attempt <- 1
  while (attempt <= max_attempts) {
    # Transform tree based on target K
    if (target_K == 0) {
      star_tree <- ape::stree(ape::Ntip(tree), type = "star")
      star_tree$tip.label <- tree$tip.label
      transformed_tree <- star_tree
    } else if (target_K == 1) {
      transformed_tree <- tree
    } else if (target_K > 1) {
      transformed_tree <- tree
      transformed_tree$edge.length <- transformed_tree$edge.length * target_K
    } else {
      # For 0 < target_K < 1
      transformed_tree <- tree
      star_tree <- ape::stree(ape::Ntip(tree), type = "star")
      star_tree$tip.label <- tree$tip.label
      mix_ratio <- target_K
      transformed_tree$edge.length <- tree$edge.length * mix_ratio +
        rep(mean(tree$edge.length), length(tree$edge.length)) * (1 - mix_ratio)
    }
    # Simulate trait data
    trait <- phytools::fastBM(transformed_tree)
    # Test phylogenetic signal
    result <- phytools::phylosig(tree, trait, method = "K", test = TRUE)
    estimated_K <- result$K
    # Check if within tolerance
    if (abs(estimated_K - target_K) <= tolerance) {
      trait_df <- data.frame(
        trait = as.numeric(trait),
        row.names = names(trait)
      )
      return(trait_df)
    }
    attempt <- attempt + 1
  }
  # If we reach here, we failed to achieve target K
  warning("Failed to achieve target K within tolerance after ",
          max_attempts, " attempts")
  return(NULL)
}
