# Test functions to generate the data

# test function_value_onemarker() ---------------------------------------------------------------------------

# Create data.frame to compare
df_one_marker  <- data.frame("count" = c(13, 4, 10, 18, 10, 6),
                             "group_id" = as.factor(rep(c("A", "B"), times = 3)),
                             "donor_id" = rep(c("donor_1", "donor_2", "donor_3"), each = 2),
                             "sample_id" = paste0("sample_", c(1,4,2,5,3,6)),
                             "marker_name" = "m1")

# Test data generation for one marker
test_that("Data generated for one marker", {
  # Set seed
  set.seed(123)
  # Generate data with function
  df_to_test <- function_value_onemarker("marker_name" = "m1",
                                         "mu0" = 10,
                                         "dispersion" = 2,
                                         "subject_effect" = 0.1,
                                         "nb_donor" = 3,
                                         "nb_cell_per_sample" = 1,
                                         "rho" = 1)
  # Test
  expect_equal(df_one_marker, df_to_test)
})

# test function_create_mock_dataset_withmarkerinfo() ---------------------------------------------------------------------------

# Create list to compare
df_exp_info <- data.frame("group_id" = as.factor(rep(c("A", "B"), times = 3)),
                          "donor_id" = rep(c("donor_1", "donor_2", "donor_3"), each = 2),
                          "sample_id" = paste0("sample_", c(1,4,2,5,3,6)))
df_raw_data <- data.frame(df_exp_info,
                          "m1" = c(14, 4, 9, 17, 10, 9),
                          "m2" = c(3, 4, 5, 14, 3, 10))
trans_data <- asinh(df_raw_data[, c("m1", "m2")]/5)
df_trans_data <- data.frame(df_exp_info,
                            trans_data)
ls_multi_marker <- list("df_info" = tibble::as_tibble(df_exp_info),
                        "DEmarkers_names" = "m2",
                        "raw_data" = tibble::as_tibble(df_raw_data),
                        "data" = tibble::as_tibble(df_trans_data))
# Create list of parameters combination for 2 markers
variation <- list(data.frame("marker_name" = "m1",
                             "nb_donor" = 3,
                             "rho" = 1,
                             "subject_effect" = 0.1,
                             "mu0" = 10,
                             "dispersion" = 3,
                             "nb_cell_per_sample" = 1),
                  data.frame("marker_name" = "m2",
                             "nb_donor" = 3,
                             "rho" = 1.1,
                             "subject_effect" = 0.1,
                             "mu0" = 10,
                             "dispersion" = 3,
                             "nb_cell_per_sample" = 1))

# Test data generated for multiple markers
test_that("Data generated for multiple marker", {
  # Set seed
  set.seed(123)
  # Generate data with function
  ls_to_test <- function_create_mock_dataset_withmarkerinfo(variation)
  # Test
  expect_equal(ls_to_test, ls_multi_marker)
})

# test function_create_mock_dataset_withmarkerinfo() ---------------------------------------------------------------------------

# Create list to compare
# Reusing data generated to test function_create_mock_dataset_withmarkerinfo()
ls_multi_sim <- list("variation" = bind_rows(variation),
                     "df_info" = tibble::as_tibble(df_exp_info),
                     "DE_markers_names" = "m2",
                     "ls_mock_data" = tibble::as_tibble(df_trans_data),
                     "ls_mock_raw_data" = tibble::as_tibble(df_raw_data))

# Test data generated for multiple markers, for multiple simulation
test_that("Data generated for multiple simulations", {
  # Set set
  set.seed(123)
  # Generate data with function
  ls_to_test <- function_apply_onesimulation_withmarkerinfo(variation)
  # Test
  expect_equal(ls_to_test, ls_multi_sim)
})
