library(testthat)

test_that("merge_dataset function tests", {

  # Prepare test data
  main_data <- data.frame(
    species = c("Abies alba", "Coussapoa trinervia", "Crataegus monogyna"),
    genus = c("Abies", "Coussapoa", "Crataegus"),
    family = c("Pinaceae", "Urticaceae", "Rosaceae"),
    LA = c(NA, 2050.24, 449.15),
    LeafN = c(13.10, 14.52, 17.46),
    Seedmass = c(53.64, NA, 95.92),
    stringsAsFactors = FALSE
  )

  additional_data <- data.frame(
    species = c("Abies alba", "Corydalis solida"),
    genus = c("Abies", "Corydalis"),
    family = c("Pinaceae", "Papaveraceae"),
    LA = c(25.58, NA),
    LMA = c(0.19, 0.2),
    PlantHeight = c(53.66, 0.14),
    stringsAsFactors = FALSE
  )

  # Test 1: Basic functionality - default priority
  result_default <- merge_dataset(main_data, additional_data)

  expect_s3_class(result_default, "data.frame")
  expect_equal(nrow(result_default), 4) # 4 unique species
  expect_true("species" %in% names(result_default))
  expect_equal(result_default$species[1], "Abies alba")

  # Test 2: Column name check
  expected_columns <- c("species", "family", "genus", "LA", "LeafN", "LMA", "PlantHeight", "Seedmass")
  expect_setequal(names(result_default), expected_columns)

  # Test 3: species column should be first
  expect_equal(names(result_default)[1], "species")

  # Test 4: Test priority "main"
  result_main <- merge_dataset(main_data, additional_data, priority = "main")
  abies_row_main <- result_main[result_main$species == "Abies alba", ]
  expect_equal(abies_row_main$family, "Pinaceae") # Both datasets have same value, use main
  expect_equal(abies_row_main$LA, 25.58) # main is NA, use additional
  expect_equal(abies_row_main$LMA, 0.19) # only additional has value

  # Test 5: Test priority "additional"
  result_additional <- merge_dataset(main_data, additional_data, priority = "additional")
  abies_row_add <- result_additional[result_additional$species == "Abies alba", ]
  expect_equal(abies_row_add$family, "Pinaceae") # Both datasets have same value
  expect_equal(abies_row_add$LA, 25.58) # main is NA, use additional

  # Test 6: Test priority "mean" for numeric data
  # Create test data with conflicting numeric values
  main_conflict <- data.frame(
    species = c("Test species"),
    value1 = c(10),
    value2 = c("text"),
    stringsAsFactors = FALSE
  )

  additional_conflict <- data.frame(
    species = c("Test species"),
    value1 = c(20),
    value2 = c("other_text"),
    stringsAsFactors = FALSE
  )

  result_mean <- merge_dataset(main_conflict, additional_conflict, priority = "mean")
  expect_equal(result_mean$value1, 15) # (10+20)/2 = 15
  expect_equal(result_mean$value2, "text") # non-numeric uses main value

  # Test 7: Missing value handling
  corydalis_row <- result_default[result_default$species == "Corydalis solida", ]
  expect_equal(corydalis_row$genus, "Corydalis")
  expect_true(is.na(corydalis_row$LA))
  expect_true(is.na(corydalis_row$LeafN))

  # Test 8: Error input validation
  expect_error(merge_dataset("not_a_dataframe", additional_data),
               "main_data and additional_data must be data frames.")

  expect_error(merge_dataset(main_data, "not_a_dataframe"),
               "main_data and additional_data must be data frames.")

  # Test 9: Missing species column error
  no_species_data <- data.frame(genus = c("Abies"), family = c("Pinaceae"))
  expect_error(merge_dataset(no_species_data, additional_data),
               "Both datasets must contain the 'species' column.")

  # Test 10: Invalid priority parameter
  expect_error(merge_dataset(main_data, additional_data, priority = "invalid"),
               "The priority parameter must be one of 'main', 'additional', or 'mean'.")

  # Test 11: Duplicate species handling - only first occurrence used
  duplicate_main <- rbind(main_data,
                          data.frame(species = "Abies alba", genus = "Different",
                                     family = "Different", LA = 999, LeafN = 999,
                                     Seedmass = 999, stringsAsFactors = FALSE))

  result_duplicate <- merge_dataset(duplicate_main, additional_data)
  abies_duplicate <- result_duplicate[result_duplicate$species == "Abies alba", ]
  expect_equal(abies_duplicate$genus, "Abies") # Should use first occurrence value
  expect_equal(abies_duplicate$LeafN, 13.10) # Should use first occurrence value

  # Test 12: Empty data frame handling
  empty_main <- data.frame(species = character(0), stringsAsFactors = FALSE)
  empty_additional <- data.frame(species = character(0), stringsAsFactors = FALSE)
  result_empty <- merge_dataset(empty_main, empty_additional)
  expect_equal(nrow(result_empty), 0)
  expect_true("species" %in% names(result_empty))

  # Test 13: Only one dataset has data
  result_only_main <- merge_dataset(main_data, empty_additional)
  expect_equal(nrow(result_only_main), 3)
  expect_setequal(result_only_main$species, main_data$species)

  # Test 14: Data type preservation
  expect_type(result_default$LA, "double")
  expect_type(result_default$species, "character")
  expect_type(result_default$family, "character")

})

# Run additional edge case tests
test_that("merge_dataset edge case tests", {

  # Test 15: Single row data
  single_main <- data.frame(species = "Single species", value = 1, stringsAsFactors = FALSE)
  single_additional <- data.frame(species = "Another species", value = 2, stringsAsFactors = FALSE)

  result_single <- merge_dataset(single_main, single_additional)
  expect_equal(nrow(result_single), 2)
  expect_setequal(result_single$species, c("Single species", "Another species"))

  # Test 16: All values are NA
  na_main <- data.frame(species = "Test", value = NA, stringsAsFactors = FALSE)
  na_additional <- data.frame(species = "Test", value = NA, stringsAsFactors = FALSE)

  result_na <- merge_dataset(na_main, na_additional)
  expect_true(is.na(result_na$value))

  # Test 17: Special characters in species names
  special_main <- data.frame(species = "Species with spaces & symbols!",
                             value = 1, stringsAsFactors = FALSE)
  special_additional <- data.frame(species = "Species with spaces & symbols!",
                                   value = 2, stringsAsFactors = FALSE)

  result_special <- merge_dataset(special_main, special_additional, priority = "mean")
  expect_equal(result_special$value, 1.5)

})

