
# Load required libraries for testing
library(testthat)
library(ape)
library(phylobase)
library(phylosignal)

# Test data setup
test_that("pnc function works correctly", {

  # Create mock phylogenetic tree
  set.seed(123)
  test_tree <- rtree(10)
  test_tree$tip.label <- paste0("sp", 1:10)

  # Create mock trait data matching tree species
  test_traits <- data.frame(
    trait1 = rnorm(10, mean = 5, sd = 2),
    trait2 = rnorm(10, mean = 10, sd = 3),
    trait3 = c(rnorm(8, mean = 2, sd = 1), NA, NA), # Include missing values
    row.names = paste0("sp", 1:10)
  )

  # Test 1: Basic functionality with default parameters
  result_basic <- pnc(test_traits, test_tree)

  expect_s3_class(result_basic, "data.frame")
  expect_true(all(c("trait", "coverage", "n_sp", "signal", "p", "significance", "method") %in% colnames(result_basic)))
  expect_equal(attr(result_basic, "methods"), "Lambda")
  expect_equal(attr(result_basic, "pca_axes"), c("PC1", "PC2"))

  # Test 2: Multiple methods specified
  result_multi <- pnc(test_traits, test_tree, methods = c("Lambda", "K", "I"))

  expect_true(nrow(result_multi) > nrow(result_basic))
  expect_true(all(c("Lambda", "K", "I") %in% result_multi$method))
  expect_equal(attr(result_multi, "methods"), c("Lambda", "K", "I"))

  # Test 3: No PCA analysis
  result_no_pca <- pnc(test_traits, test_tree, pca_axes = NULL)

  expect_false(any(grepl("PC", result_no_pca$trait)))
  expect_null(attr(result_no_pca, "pca_results"))

  # Test 4: Custom significance levels
  custom_sig <- c(0.001, 0.01, 0.05, 0.1)
  result_custom_sig <- pnc(test_traits, test_tree, sig_levels = custom_sig)

  expect_equal(attr(result_custom_sig, "sig_levels"), custom_sig)

  # Test 5: Custom number of repetitions
  result_custom_reps <- pnc(test_traits, test_tree, reps = 99)

  expect_equal(attr(result_custom_reps, "reps"), 99)

  # Test 6: Verbose mode off
  expect_silent(pnc(test_traits, test_tree, verbose = FALSE))
})

# Test error handling and edge cases
test_that("pnc handles edge cases and errors correctly", {

  # Setup test data
  set.seed(123)
  test_tree <- rtree(10)
  test_tree$tip.label <- paste0("sp", 1:10)

  test_traits <- data.frame(
    trait1 = rnorm(10),
    trait2 = rnorm(10),
    row.names = paste0("sp", 1:10)
  )

  # Test invalid method specification
  expect_error(pnc(test_traits, test_tree, methods = "invalid_method"))

  # Test with insufficient species (< 4)
  small_tree <- rtree(3)
  small_tree$tip.label <- paste0("sp", 1:3)
  small_traits <- data.frame(
    trait1 = rnorm(3),
    row.names = paste0("sp", 1:3)
  )

  result_small <- pnc(small_traits, small_tree, verbose = FALSE)
  expect_true(all(is.na(result_small$signal)))
  expect_true(all(is.na(result_small$p)))

  # Test with mismatched species names
  mismatch_traits <- test_traits
  rownames(mismatch_traits) <- paste0("other_sp", 1:10)

  result_mismatch <- pnc(mismatch_traits, test_tree, verbose = FALSE)
  expect_true(all(is.na(result_mismatch$signal)))

  # Test with all missing trait data
  na_traits <- data.frame(
    trait1 = rep(NA, 10),
    row.names = paste0("sp", 1:10)
  )

  result_na <- pnc(na_traits, test_tree, verbose = FALSE)
  expect_true(all(is.na(result_na$signal)))
})

# Test PCA functionality
test_that("PCA analysis works correctly", {

  # Setup test data with sufficient observations for PCA
  set.seed(123)
  test_tree <- rtree(15)
  test_tree$tip.label <- paste0("sp", 1:15)

  # Create correlated traits for meaningful PCA
  trait_base <- rnorm(15, mean = 5, sd = 2)
  test_traits <- data.frame(
    trait1 = trait_base + rnorm(15, sd = 0.5),
    trait2 = trait_base * 1.5 + rnorm(15, sd = 1),
    trait3 = trait_base * -0.8 + rnorm(15, sd = 0.8),
    trait4 = rnorm(15, mean = 10, sd = 2),
    row.names = paste0("sp", 1:15)
  )

  # Test successful PCA
  result_pca <- pnc(test_traits, test_tree, pca_axes = c("PC1", "PC2", "PC3"))

  expect_true(any(grepl("PC", result_pca$trait)))
  expect_false(attr(result_pca, "pca_failed"))
  expect_s3_class(attr(result_pca, "pca_results"), "data.frame")

  # Test PCA with insufficient complete cases
  sparse_traits <- test_traits
  sparse_traits[1:12, ] <- NA # Leave only 3 complete cases

  result_sparse <- pnc(sparse_traits, test_tree, verbose = FALSE)
  expect_true(attr(result_sparse, "pca_failed"))
})

# Test significance marking function
test_that("Significance levels are correctly assigned", {

  set.seed(123)
  test_tree <- rtree(8)
  test_tree$tip.label <- paste0("sp", 1:8)

  test_traits <- data.frame(
    trait1 = rnorm(8),
    row.names = paste0("sp", 1:8)
  )

  result <- pnc(test_traits, test_tree, sig_levels = c(0.001, 0.01, 0.05))

  # Check that significance column contains expected values
  expect_true(all(result$significance %in% c("***", "**", "*", "ns", "")))
})

# Integration test with different phylogenetic signal methods
test_that("All phylogenetic signal methods produce valid output", {

  set.seed(123)
  test_tree <- rtree(12)
  test_tree$tip.label <- paste0("sp", 1:12)

  # Create trait with some phylogenetic signal
  test_trait_values <- rnorm(12)
  names(test_trait_values) <- test_tree$tip.label

  test_traits <- data.frame(
    phylo_trait = test_trait_values,
    row.names = names(test_trait_values)
  )

  # Test all available methods
  all_methods <- c("Lambda", "K", "K.star", "I", "Cmean")
  result_all <- pnc(test_traits, test_tree, methods = all_methods, pca_axes = NULL, verbose = FALSE)

  expect_equal(length(unique(result_all$method)), length(all_methods))
  expect_true(all(all_methods %in% result_all$method))

  # Check that each method produces results
  for(method in all_methods) {
    method_results <- result_all[result_all$method == method, ]
    expect_equal(nrow(method_results), 1)
    expect_true(!is.null(method_results$signal))
    expect_true(!is.null(method_results$p))
  }
})


# Helper function for creating test datasets
create_multi_community_data <- function(n_plots = 5, n_species = 15, n_traits = 3, missing_prop = 0.1) {
  set.seed(123)

  # Create phylogenetic tree
  tree <- rtree(n_species)
  tree$tip.label <- paste0("sp", 1:n_species)

  # Create community matrix
  com_matrix <- matrix(rpois(n_plots * n_species, lambda = 2),
                       nrow = n_plots, ncol = n_species)
  colnames(com_matrix) <- tree$tip.label
  rownames(com_matrix) <- paste0("plot", sprintf("%03d", 1:n_plots))

  # Create trait matrix
  traits <- matrix(rnorm(n_species * n_traits), nrow = n_species, ncol = n_traits)
  colnames(traits) <- paste0("trait", 1:n_traits)
  rownames(traits) <- tree$tip.label

  # Introduce missing values
  if (missing_prop > 0) {
    n_missing <- floor(n_species * n_traits * missing_prop)
    missing_indices <- sample(1:(n_species * n_traits), n_missing)
    traits[missing_indices] <- NA
  }

  traits_df <- as.data.frame(traits)

  return(list(tree = tree, com = com_matrix, traits = traits_df))
}

# Performance test
test_that("compnc function performance is acceptable", {

  # Create larger dataset
  large_data <- create_multi_community_data(n_plots = 10, n_species = 50, n_traits = 5)

  # Measure execution time
  start_time <- Sys.time()
  result <- compnc(large_data$com, large_data$traits, large_data$tree, verbose = FALSE)
  end_time <- Sys.time()

  execution_time <- as.numeric(difftime(end_time, start_time, units = "secs"))

  # Should complete within reasonable time (adjust threshold as needed)
  expect_lt(execution_time, 60)  # Less than 60 seconds

  # Results should still be valid
  expect_s3_class(result, "data.frame")
  expect_true(nrow(result) > 0)
  expect_equal(attr(result, "total_plots"), 10)
})






