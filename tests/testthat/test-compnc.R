
# Load required libraries for testing
library(testthat)
library(ape)
library(phylobase)
library(phylosignal)

# Test data setup
test_that("compnc function works correctly", {

  # Create mock phylogenetic tree
  set.seed(123)
  test_tree <- rtree(15)
  test_tree$tip.label <- paste0("sp", 1:15)

  # Create mock community matrix (plots x species)
  test_com <- matrix(rpois(45, lambda = 2), nrow = 3, ncol = 15)
  colnames(test_com) <- test_tree$tip.label
  rownames(test_com) <- paste0("plot", 1:3)

  # Create mock trait data matching tree species
  test_traits <- data.frame(
    trait1 = rnorm(15, mean = 5, sd = 2),
    trait2 = rnorm(15, mean = 10, sd = 3),
    trait3 = c(rnorm(12, mean = 2, sd = 1), NA, NA, NA), # Include missing values
    row.names = test_tree$tip.label
  )

  # Test 1: Basic functionality with default parameters
  result_basic <- compnc(test_com, test_traits, test_tree)

  expect_s3_class(result_basic, "data.frame")
  expect_true(all(c("plot", "trait", "coverage", "n_sp", "signal", "p", "significance", "method", "n_sp_in_plot") %in% colnames(result_basic)))
  expect_equal(attr(result_basic, "methods"), c("Lambda", "K", "K.star", "I", "Cmean"))
  expect_equal(attr(result_basic, "pca_axes"), c("PC1", "PC2"))

  # Test 2: Multiple methods specified
  result_multi <- compnc(test_com, test_traits, test_tree, methods = c("Lambda", "K", "I"))

  expect_true(all(c("Lambda", "K", "I") %in% result_multi$method))
  expect_equal(attr(result_multi, "methods"), c("Lambda", "K", "I"))

  # Test 3: No PCA analysis
  result_no_pca <- compnc(test_com, test_traits, test_tree, pca_axes = NULL)

  expect_false(any(grepl("PC", result_no_pca$trait)))

  # Test 4: Custom significance levels
  custom_sig <- c(0.001, 0.01, 0.05, 0.1)
  result_custom_sig <- compnc(test_com, test_traits, test_tree, sig_levels = custom_sig)

  expect_equal(attr(result_custom_sig, "sig_levels"), custom_sig)

  # Test 5: Custom minimum abundance threshold
  result_min_abund <- compnc(test_com, test_traits, test_tree, min_abundance = 1)

  expect_equal(attr(result_min_abund, "min_abundance"), 1)

  # Test 6: Custom number of repetitions
  result_custom_reps <- compnc(test_com, test_traits, test_tree, reps = 99)

  expect_equal(attr(result_custom_reps, "reps"), 99)

  # Test 7: Verbose mode off
  expect_silent(compnc(test_com, test_traits, test_tree, verbose = FALSE))
})

# Test error handling and edge cases
test_that("compnc handles edge cases and errors correctly", {

  # Setup test data
  set.seed(123)
  test_tree <- rtree(10)
  test_tree$tip.label <- paste0("sp", 1:10)

  test_com <- matrix(rpois(30, lambda = 2), nrow = 3, ncol = 10)
  colnames(test_com) <- test_tree$tip.label
  rownames(test_com) <- paste0("plot", 1:3)

  test_traits <- data.frame(
    trait1 = rnorm(10),
    trait2 = rnorm(10),
    row.names = test_tree$tip.label
  )

  # Test invalid method specification
  expect_error(compnc(test_com, test_traits, test_tree, methods = "invalid_method"))

  # Test with insufficient species in community (< 4)
  small_com <- matrix(c(1, 0, 0, 2, 1, 0), nrow = 2, ncol = 3)
  colnames(small_com) <- test_tree$tip.label[1:3]
  rownames(small_com) <- paste0("plot", 1:2)

  result_small <- compnc(small_com, test_traits, test_tree, verbose = FALSE)

  # Should return empty result or with NA values due to insufficient species
  expect_s3_class(result_small, "data.frame")

  # Test with mismatched species names between community and traits
  mismatch_com <- test_com
  colnames(mismatch_com) <- paste0("other_sp", 1:10)

  result_mismatch <- compnc(mismatch_com, test_traits, test_tree, verbose = FALSE)
  expect_s3_class(result_mismatch, "data.frame")

  # Test with all missing trait data
  na_traits <- data.frame(
    trait1 = rep(NA, 10),
    row.names = test_tree$tip.label
  )

  result_na <- compnc(test_com, na_traits, test_tree, verbose = FALSE)
  expect_s3_class(result_na, "data.frame")
})

# Test community matrix handling
test_that("Community matrix processing works correctly", {

  # Setup test data
  set.seed(123)
  test_tree <- rtree(12)
  test_tree$tip.label <- paste0("sp", 1:12)

  test_traits <- data.frame(
    trait1 = rnorm(12, mean = 5, sd = 2),
    trait2 = rnorm(12, mean = 10, sd = 3),
    row.names = test_tree$tip.label
  )

  # Test 1: Community matrix without row names (auto-generated)
  test_com_no_names <- matrix(rpois(48, lambda = 2), nrow = 4, ncol = 12)
  colnames(test_com_no_names) <- test_tree$tip.label

  result_auto_names <- compnc(test_com_no_names, test_traits, test_tree, verbose = FALSE)

  expect_true(all(grepl("plot", unique(result_auto_names$plot))))
  expect_equal(attr(result_auto_names, "total_plots"), 4)

  # Test 2: Different abundance patterns
  sparse_com <- matrix(0, nrow = 3, ncol = 12)
  colnames(sparse_com) <- test_tree$tip.label
  rownames(sparse_com) <- paste0("plot", 1:3)
  # Only add species to first plot
  sparse_com[1, 1:8] <- c(5, 3, 7, 2, 4, 6, 1, 8)

  result_sparse <- compnc(sparse_com, test_traits, test_tree, verbose = FALSE)

  # Should only analyze the first plot with sufficient species
  analyzed_plots <- length(unique(result_sparse$plot[!is.na(result_sparse$signal)]))
  expect_true(analyzed_plots <= 1)

  # Test 3: Minimum abundance filtering
  result_min_abund <- compnc(sparse_com, test_traits, test_tree, min_abundance = 3, verbose = FALSE)

  # Should further filter species based on minimum abundance
  expect_s3_class(result_min_abund, "data.frame")
})

# Test PCA functionality for multiple communities
test_that("PCA analysis works correctly across communities", {

  # Setup test data with sufficient observations for PCA
  set.seed(123)
  test_tree <- rtree(20)
  test_tree$tip.label <- paste0("sp", 1:20)

  # Create community matrix with good coverage
  test_com <- matrix(rpois(100, lambda = 3), nrow = 5, ncol = 20)
  colnames(test_com) <- test_tree$tip.label
  rownames(test_com) <- paste0("plot", 1:5)

  # Create correlated traits for meaningful PCA
  trait_base <- rnorm(20, mean = 5, sd = 2)
  test_traits <- data.frame(
    trait1 = trait_base + rnorm(20, sd = 0.5),
    trait2 = trait_base * 1.5 + rnorm(20, sd = 1),
    trait3 = trait_base * -0.8 + rnorm(20, sd = 0.8),
    trait4 = rnorm(20, mean = 10, sd = 2),
    row.names = test_tree$tip.label
  )

  # Test successful PCA across multiple communities
  result_pca <- compnc(test_com, test_traits, test_tree, pca_axes = c("PC1", "PC2", "PC3"), verbose = FALSE)

  expect_true(any(grepl("PC", result_pca$trait)))
  expect_equal(attr(result_pca, "pca_axes"), c("PC1", "PC2", "PC3"))
})

# Test significance marking function
test_that("Significance levels are correctly assigned", {

  set.seed(123)
  test_tree <- rtree(15)
  test_tree$tip.label <- paste0("sp", 1:15)

  test_com <- matrix(rpois(45, lambda = 3), nrow = 3, ncol = 15)
  colnames(test_com) <- test_tree$tip.label
  rownames(test_com) <- paste0("plot", 1:3)

  test_traits <- data.frame(
    trait1 = rnorm(15),
    row.names = test_tree$tip.label
  )

  result <- compnc(test_com, test_traits, test_tree,
                   sig_levels = c(0.001, 0.01, 0.05),
                   pca_axes = NULL,
                   verbose = FALSE)

  # Check that significance column contains expected values
  expect_true(all(result$significance %in% c("***", "**", "*", "ns", "")))

  # Test coverage calculation format
  expect_true(all(grepl("%", result$coverage[!is.na(result$coverage)])))
})

# Integration test with all phylogenetic signal methods
test_that("All phylogenetic signal methods produce valid output", {

  set.seed(123)
  test_tree <- rtree(18)
  test_tree$tip.label <- paste0("sp", 1:18)

  # Create community matrix
  test_com <- matrix(rpois(72, lambda = 2), nrow = 4, ncol = 18)
  colnames(test_com) <- test_tree$tip.label
  rownames(test_com) <- paste0("plot", 1:4)

  # Create trait with some phylogenetic signal
  test_traits <- data.frame(
    phylo_trait = rnorm(18),
    row.names = test_tree$tip.label
  )

  # Test all available methods
  all_methods <- c("Lambda", "K", "K.star", "I", "Cmean")
  result_all <- compnc(test_com, test_traits, test_tree,
                       methods = all_methods,
                       pca_axes = NULL,
                       verbose = FALSE)

  expect_true(all(all_methods %in% result_all$method))

  # Check that each method produces results for valid communities
  for(method in all_methods) {
    method_results <- result_all[result_all$method == method, ]
    expect_true(nrow(method_results) > 0)
  }

  # Check that n_sp_in_plot is correctly calculated
  expect_true(all(result_all$n_sp_in_plot >= result_all$n_sp, na.rm = TRUE))
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







