library(testthat)

test_that("extract_traits function tests", {

  # Prepare test dataset
  test_dataset <- data.frame(
    species = c("Acaena novae-zelandiae", "Acaena novae-zelandiae", "Adiantum capillus-veneris",
                "Zuelania guidonia", "Betula pendula", "Betula pubescens"),
    genus = c("Acaena", "Acaena", "Adiantum", "Zuelania", "Betula", "Betula"),
    family = c("Rosaceae", "Rosaceae", "Pteridaceae", "Salicaceae", "Betulaceae", "Betulaceae"),
    LA = c(25.5, 30.0, 150.2, NA, 45.8, 52.1),
    LMA = c(0.15, 0.18, 0.12, 0.25, NA, 0.22),
    LeafN = c(12.5, 13.2, 8.9, 15.6, 11.8, 12.4),
    PlantHeight = c(0.3, 0.35, 0.4, 25.5, 18.2, 20.1),
    SeedMass = c(2.1, 2.3, NA, 125.6, 0.8, 0.9),
    NonNumeric = c("type1", "type2", "type3", "type4", "type5", "type6"),
    stringsAsFactors = FALSE
  )

  # Test 1: Basic functionality - species level, all traits
  species_list <- c("Acaena novae-zelandiae", "Adiantum capillus-veneris")
  result_species_all <- extract_traits(species_list, test_dataset, rank = "species")

  expect_s3_class(result_species_all, "data.frame")
  expect_equal(nrow(result_species_all), 2)
  expect_equal(rownames(result_species_all), species_list)
  expect_true("LA" %in% colnames(result_species_all))
  expect_true("NonNumeric" %in% colnames(result_species_all))

  # Test 2: Species level with specific traits
  specific_traits <- c("LA", "LMA", "LeafN")
  result_species_specific <- extract_traits(species_list, test_dataset,
                                            rank = "species", traits = specific_traits)

  expect_equal(ncol(result_species_specific), 3)
  expect_setequal(colnames(result_species_specific), specific_traits)
  expect_equal(result_species_specific["Acaena novae-zelandiae", "LA"], 25.5)

  # Test 3: Genus level - mean calculation for numeric traits
  genus_list <- c("Acaena", "Betula")
  result_genus <- extract_traits(genus_list, test_dataset, rank = "genus",
                                 traits = c("LA", "LMA", "NonNumeric"))

  expect_equal(nrow(result_genus), 2)
  expect_equal(result_genus["Acaena", "LA"], mean(c(25.5, 30.0), na.rm = TRUE))
  expect_equal(result_genus["Betula", "LMA"], mean(c(0.22), na.rm = TRUE)) # NA excluded
  expect_equal(result_genus["Acaena", "NonNumeric"], "type1") # First occurrence for non-numeric

  # Test 4: Family level - mean calculation
  family_list <- c("Rosaceae", "Betulaceae")
  result_family <- extract_traits(family_list, test_dataset, rank = "family",
                                  traits = c("LA", "PlantHeight"))

  expect_equal(nrow(result_family), 2)
  expect_equal(result_family["Rosaceae", "LA"], mean(c(25.5, 30.0), na.rm = TRUE))
  expect_equal(result_family["Betulaceae", "PlantHeight"], mean(c(18.2, 20.1), na.rm = TRUE))

  # Test 5: Error handling - invalid rank
  expect_error(extract_traits(species_list, test_dataset, rank = "invalid"),
               "rank must be one of 'species', 'genus', or 'family'")

  # Test 6: Error handling - no valid traits specified
  expect_error(extract_traits(species_list, test_dataset,
                              traits = c("invalid_trait1", "invalid_trait2")),
               "No valid traits specified")

  # Test 7: Warning for invalid traits (capture output)
  expect_output(extract_traits(species_list, test_dataset,
                               traits = c("LA", "invalid_trait")),
                "Warning: The following traits are not available in the dataset:")

  # Test 8: Missing taxa handling
  missing_species <- c("Acaena novae-zelandiae", "Missing species")
  expect_output(extract_traits(missing_species, test_dataset, rank = "species"),
                "The missing species in the dataset:")

  result_missing <- extract_traits(missing_species, test_dataset, rank = "species")
  expect_true(all(is.na(result_missing["Missing species", ])))

  # Test 9: Empty species list
  empty_result <- extract_traits(character(0), test_dataset, rank = "species")
  expect_equal(nrow(empty_result), 0)
  expect_true(ncol(empty_result) > 0) # Should have trait columns

  # Test 10: Single species
  single_result <- extract_traits("Acaena novae-zelandiae", test_dataset, rank = "species")
  expect_equal(nrow(single_result), 1)
  expect_equal(rownames(single_result), "Acaena novae-zelandiae")

  # Test 11: All NA values for a trait
  na_dataset <- test_dataset
  na_dataset$LA[na_dataset$genus == "Acaena"] <- NA
  result_na <- extract_traits("Acaena", na_dataset, rank = "genus", traits = "LA")
  expect_true(is.na(result_na["Acaena", "LA"]))

  # Test 12: Mixed data types handling
  mixed_result <- extract_traits("Acaena", test_dataset, rank = "genus",
                                 traits = c("LA", "NonNumeric"))
  expect_type(mixed_result$LA, "double")
  expect_type(mixed_result$NonNumeric, "character")

  # Test 13: Duplicate species handling (first occurrence)
  dup_result <- extract_traits("Acaena novae-zelandiae", test_dataset, rank = "species")
  expect_equal(dup_result["Acaena novae-zelandiae", "LA"], 25.5) # First occurrence
  expect_equal(dup_result["Acaena novae-zelandiae", "NonNumeric"], "type1")

  # Test 14: NULL traits parameter (all traits)
  all_traits_result <- extract_traits("Acaena novae-zelandiae", test_dataset,
                                      rank = "species", traits = NULL)
  expected_traits <- setdiff(colnames(test_dataset), c("species", "genus", "family"))
  expect_setequal(colnames(all_traits_result), expected_traits)

})

test_that("extract_traits edge cases", {

  # Prepare minimal test dataset
  minimal_dataset <- data.frame(
    species = c("Species1", "Species2"),
    genus = c("Genus1", "Genus2"),
    family = c("Family1", "Family2"),
    trait1 = c(10, 20),
    stringsAsFactors = FALSE
  )

  # Test 15: Dataset with minimum required columns
  result_minimal <- extract_traits("Species1", minimal_dataset, rank = "species")
  expect_equal(ncol(result_minimal), 1)
  expect_equal(colnames(result_minimal), "trait1")
  expect_equal(result_minimal["Species1", "trait1"], 10)

  # Test 16: Dataset with only taxonomic columns (no traits)
  no_traits_dataset <- data.frame(
    species = c("Species1"),
    genus = c("Genus1"),
    family = c("Family1"),
    stringsAsFactors = FALSE
  )

  result_no_traits <- extract_traits("Species1", no_traits_dataset, rank = "species")
  expect_equal(ncol(result_no_traits), 0)
  expect_equal(nrow(result_no_traits), 1)

  # Test 17: Special characters in taxa names
  special_dataset <- data.frame(
    species = c("Species with spaces & symbols!"),
    genus = c("Genus-name"),
    family = c("Family_name"),
    trait1 = c(5.5),
    stringsAsFactors = FALSE
  )

  result_special <- extract_traits("Species with spaces & symbols!", special_dataset,
                                   rank = "species")
  expect_equal(result_special["Species with spaces & symbols!", "trait1"], 5.5)

  # Test 18: Case sensitivity
  case_dataset <- data.frame(
    species = c("species lower", "SPECIES UPPER"),
    genus = c("genus lower", "GENUS UPPER"),
    family = c("family lower", "FAMILY UPPER"),
    trait1 = c(1, 2),
    stringsAsFactors = FALSE
  )

  # Should not match due to case difference
  expect_output(extract_traits("Species Lower", case_dataset, rank = "species"),
                "The missing species in the dataset:")

  # Test 19: NaN handling
  nan_dataset <- data.frame(
    species = c("TestSpecies"),
    genus = c("TestGenus"),
    family = c("TestFamily"),
    trait1 = c(NaN),
    stringsAsFactors = FALSE
  )

  result_nan <- extract_traits("TestSpecies", nan_dataset, rank = "species")
  expect_true(is.na(result_nan["TestSpecies", "trait1"]))

  # Test 20: Large dataset performance (basic check)
  large_dataset <- data.frame(
    species = paste0("Species", 1:1000),
    genus = paste0("Genus", rep(1:100, each = 10)),
    family = paste0("Family", rep(1:10, each = 100)),
    trait1 = rnorm(1000),
    trait2 = runif(1000),
    stringsAsFactors = FALSE
  )

  large_result <- extract_traits(c("Species1", "Species500", "Species1000"),
                                 large_dataset, rank = "species")
  expect_equal(nrow(large_result), 3)
  expect_equal(ncol(large_result), 2)

})
