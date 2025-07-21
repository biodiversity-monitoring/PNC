library(testthat)

test_that("coverage function tests", {

  # Prepare the test dataset
  test_data <- data.frame(
    PlantHeight = c(1.2, 1.5, NA, 2.1, 1.8),
    LDMC = c(0.5, NA, 0.8, 1.2, 0.9),
    LA = c(15.2, 18.5, 12.3, NA, 16.7),
    SeedMass = c(NA, NA, NA, NA, NA)
  )

  # Test 1: Basic Function Test
  expect_output(result <- coverage(test_data), "Overall trait coverage rate:")

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 5)
  expect_equal(ncol(result), 4)
  expect_setequal(colnames(result), c("Trait", "Available_count", "Missing_count", "Trait_coverage_rate"))

  # Test 2: Column Name Check
  expected_traits <- c("PlantHeight", "LDMC", "LA", "SeedMass", "All")
  expect_setequal(result$Trait, expected_traits)

  # Test 3: Specific numerical verification
  plant_height_row <- result[result$Trait == "PlantHeight", ]
  expect_equal(plant_height_row$Available_count, 4)
  expect_equal(plant_height_row$Missing_count, 1)
  expect_equal(plant_height_row$Trait_coverage_rate, "80 %")

  ldmc_row <- result[result$Trait == "LDMC", ]
  expect_equal(ldmc_row$Available_count, 4)
  expect_equal(ldmc_row$Missing_count, 1)
  expect_equal(ldmc_row$Trait_coverage_rate, "80 %")

  # Test 4: All missing traits
  seed_mass_row <- result[result$Trait == "SeedMass", ]
  expect_equal(seed_mass_row$Available_count, 0)
  expect_equal(seed_mass_row$Missing_count, 5)
  expect_equal(seed_mass_row$Trait_coverage_rate, "0 %")

  # Test 5: "All" line (complete case) verification
  all_row <- result[result$Trait == "All", ]
  complete_cases_count <- sum(complete.cases(test_data))
  expect_equal(all_row$Available_count, complete_cases_count)
  expect_equal(all_row$Missing_count, 5 - complete_cases_count)

  # Test 6: Percentage Format Check
  expect_true(all(grepl(" %$", result$Trait_coverage_rate)))

  # Test 7: Numeric Type Checking
  expect_type(result$Available_count, "integer")
  expect_type(result$Missing_count, "integer")
  expect_type(result$Trait_coverage_rate, "character")
  expect_type(result$Trait, "character")
})

test_that("coverage function edge cases", {

  # Test 8: No missing data
  complete_data <- data.frame(
    trait1 = c(1, 2, 3),
    trait2 = c(4, 5, 6),
    trait3 = c(7, 8, 9)
  )

  expect_output(result_complete <- coverage(complete_data), "Overall trait coverage rate: 100 %")

  expect_equal(nrow(result_complete), 4)
  expect_true(all(result_complete$Trait_coverage_rate[1:3] == "100 %"))
  expect_equal(result_complete$Available_count[result_complete$Trait == "All"], 3)
  expect_equal(result_complete$Missing_count[result_complete$Trait == "All"], 0)

  # Test 9: All missing data
  all_na_data <- data.frame(
    trait1 = c(NA, NA, NA),
    trait2 = c(NA, NA, NA)
  )

  expect_output(result_na <- coverage(all_na_data), "Overall trait coverage rate: 0 %")

  expect_true(all(result_na$Trait_coverage_rate[1:2] == "0 %"))
  expect_equal(result_na$Available_count[result_na$Trait == "All"], 0)
  expect_equal(result_na$Missing_count[result_na$Trait == "All"], 3)

  # Test 10: Single-column data
  single_col_data <- data.frame(
    only_trait = c(1, NA, 3, NA, 5)
  )

  expect_output(result_single <- coverage(single_col_data), "Overall trait coverage rate: 60 %")

  expect_equal(nrow(result_single), 2)  # 1个性状 + 1个"All"行
  expect_equal(result_single$Available_count[1], 3)
  expect_equal(result_single$Missing_count[1], 2)
  expect_equal(result_single$Trait_coverage_rate[1], "60 %")

  # Test 11: Single-row data
  single_row_data <- data.frame(
    trait1 = 1,
    trait2 = NA,
    trait3 = 3
  )

  expect_output(result_single_row <- coverage(single_row_data), "Overall trait coverage rate:")

  expect_equal(nrow(result_single_row), 4)
  expect_equal(result_single_row$Available_count[result_single_row$Trait == "trait1"], 1)
  expect_equal(result_single_row$Missing_count[result_single_row$Trait == "trait2"], 1)

  # Test 12: Empty Data Frame (Should generate error or special handling)
  empty_data <- data.frame()
  expect_error(coverage(empty_data))

  # Test 13: Complete data with only one observation
  minimal_complete <- data.frame(
    trait1 = 1,
    trait2 = 2
  )

  expect_output(result_minimal <- coverage(minimal_complete), "Overall trait coverage rate: 100 %")
  expect_equal(result_minimal$Available_count[result_minimal$Trait == "All"], 1)

  # Test 14: Data containing special values
  special_data <- data.frame(
    trait1 = c(0, -1, Inf, -Inf, 1),
    trait2 = c(NA, 0, -0, 1e-10, 1e10)
  )

  expect_output(result_special <- coverage(special_data))
  expect_equal(result_special$Available_count[result_special$Trait == "trait1"], 5)
  expect_equal(result_special$Missing_count[result_special$Trait == "trait2"], 1)
})

test_that("coverage function mathematical accuracy", {

  # Test 15: Precise mathematical calculation verification
  precise_data <- data.frame(
    trait1 = c(1, 2, NA),      # 2/3 = 66.67%
    trait2 = c(NA, NA, 3),     # 1/3 = 33.33%
    trait3 = c(1, 2, 3)        # 3/3 = 100%
  )

  expect_output(result_precise <- coverage(precise_data))

  expect_equal(result_precise$Trait_coverage_rate[result_precise$Trait == "trait1"], "66.67 %")
  expect_equal(result_precise$Trait_coverage_rate[result_precise$Trait == "trait2"], "33.33 %")
  expect_equal(result_precise$Trait_coverage_rate[result_precise$Trait == "trait3"], "100 %")

  expect_equal(result_precise$Available_count[result_precise$Trait == "All"], 0)
  expect_equal(result_precise$Trait_coverage_rate[result_precise$Trait == "All"], "0 %")

  # Test 16: Verification of Overall coverage rate calculation
  expect_output(coverage(precise_data), "Overall trait coverage rate: 66.67 %")

  # Test 17: Data with complete cases
  complete_case_data <- data.frame(
    trait1 = c(1, NA, 3),
    trait2 = c(2, 4, NA),
    trait3 = c(5, 6, 7)
  )

  expect_output(result_complete_case <- coverage(complete_case_data))
  expect_equal(result_complete_case$Available_count[result_complete_case$Trait == "All"], 1)
  expect_equal(result_complete_case$Trait_coverage_rate[result_complete_case$Trait == "All"], "33.33 %")
})

test_that("coverage function data types", {

  # Test 18: Mixed data types
  mixed_data <- data.frame(
    numeric_trait = c(1.5, 2.3, NA),
    integer_trait = c(1L, 2L, 3L),
    character_trait = c("a", NA, "c"),
    logical_trait = c(TRUE, FALSE, NA),
    factor_trait = factor(c("level1", "level2", NA)),
    stringsAsFactors = FALSE
  )

  expect_output(result_mixed <- coverage(mixed_data))

  expect_equal(nrow(result_mixed), 6)  # 5个性状 + 1个"All"行
  expect_equal(result_mixed$Available_count[result_mixed$Trait == "numeric_trait"], 2)
  expect_equal(result_mixed$Available_count[result_mixed$Trait == "character_trait"], 2)
  expect_equal(result_mixed$Available_count[result_mixed$Trait == "logical_trait"], 2)
  expect_equal(result_mixed$Available_count[result_mixed$Trait == "factor_trait"], 2)

  # Test 19: Basic Functions of Extremely Large Datasets
  large_data <- data.frame(
    trait1 = sample(c(1:100, NA), 1000, replace = TRUE),
    trait2 = sample(c(1:100, NA), 1000, replace = TRUE)
  )

  expect_output(result_large <- coverage(large_data))
  expect_equal(nrow(result_large), 3)
  expect_true(all(result_large$Available_count + result_large$Missing_count == c(1000, 1000, 1000)))
})

test_that("coverage function output format", {

  # Test 20: Consistency of Output format
  format_data <- data.frame(
    trait1 = c(1, 2, 3, 4, 5),
    trait2 = c(NA, 2, 3, 4, 5)
  )

  expect_output(result_format <- coverage(format_data), "Overall trait coverage rate:")

  percentages <- gsub(" %", "", result_format$Trait_coverage_rate)
  expect_true(all(grepl("^[0-9]+(\\.[0-9]{1,2})?$", percentages)))

  expect_true(is.character(rownames(result_format)) || is.null(rownames(result_format)))

  # Test 21: Console output capture
  expect_output(coverage(format_data), "Overall trait coverage rate: 90 %")
})

