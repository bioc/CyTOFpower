# Test functions that run models for one simulation

# Create list of parameters combination for 3 markers
variation <- list(data.frame("marker_name" = "m1",
                             "nb_donor" = 3,
                             "rho" = 1,
                             "subject_effect" = 0.1,
                             "mu0" = 10,
                             "dispersion" = 3,
                             "nb_cell_per_sample" = 3),
                  data.frame("marker_name" = "m2",
                             "nb_donor" = 3,
                             "rho" = 1.1,
                             "subject_effect" = 0.1,
                             "mu0" = 10,
                             "dispersion" = 3,
                             "nb_cell_per_sample" = 3),
                  data.frame("marker_name" = "m3",
                             "nb_donor" = 3,
                             "rho" = 1,
                             "subject_effect" = 0.1,
                             "mu0" = 10,
                             "dispersion" = 3,
                             "nb_cell_per_sample" = 3))
# Create a mock dataset
df_exp_info <- data.frame("group_id" = as.factor(rep(c("A", "B"), times = 3)),
                          "donor_id" = rep(c("donor_1", "donor_2", "donor_3"), each = 2),
                          "sample_id" = paste0("sample_", c(1,4,2,5,3,6)))
df_exp_info_rep <- df_exp_info[rep(seq_len(nrow(df_exp_info)), each = 3), ]
df_raw_data <- data.frame(df_exp_info,
                          "m1" = c(14, 4, 9, 17, 10, 9, 4, 2, 5, 2, 14, 1, 11, 4, 10, 9, 22, 16),
                          "m2" = c(10, 5, 16, 12, 12, 13, 10, 3, 11, 6, 5, 6, 6, 5, 16, 8, 17, 8),
                          "m3" = c(8, 5, 6, 3, 6, 10, 7, 2, 12, 7, 4, 14, 12, 24, 5, 12, 20, 20))
trans_data <- asinh(df_raw_data[, c("m1", "m2", "m3")]/5)
df_trans_data <- data.frame(df_exp_info,
                            trans_data)
ls_3markers <- list("variation" = bind_rows(variation),
                    "df_info" = tibble::as_tibble(df_exp_info),
                    "DE_markers_names" = "m2",
                    "ls_mock_data" = tibble::as_tibble(df_trans_data),
                    "ls_mock_raw_data" = tibble::as_tibble(df_raw_data))

# Test function_to_compute_model_computation_onesimulation_modelchoice() ---------------------------------------------------------------------------

# Test
test_that("Run models", {
  # Set seed
  set.seed(123)
  # Run model - cytoglmm
  res_cytoglmm <- function_to_compute_model_computation_onesimulation_modelchoice(
    onevariation = ls_3markers,
    model = "cytoglmm")

  # Test output
  # Is it a list?
  expect_type(res_cytoglmm, "list")
  # Does it have the right structure?
  expect_identical(names(res_cytoglmm),
                   c("model", "marker_id", "p_val", "p_adj", "truth"))

  # Run model - testDS_limma_random
  res_limma_random <- function_to_compute_model_computation_onesimulation_modelchoice(
    onevariation = ls_3markers,
    model = "testDS_limma_random")
  # Test output
  # Is it a list?
  expect_type(res_limma_random, "list")
  # Does it have the right structure?
  expect_identical(names(res_limma_random),
                   c("model", "marker_id", "p_val", "p_adj", "truth"))

  # Run model - testDS_limma_fixed
  res_limma_fixed <- function_to_compute_model_computation_onesimulation_modelchoice(
    onevariation = ls_3markers,
    model = "testDS_limma_fixed")
  # Test output
  # Is it a list?
  expect_type(res_limma_fixed, "list")
  # Does it have the right structure?
  expect_identical(names(res_limma_fixed),
                   c("model", "marker_id", "p_val", "p_adj", "truth"))

  # Run model - testDS_lmm
  res_lmm <- function_to_compute_model_computation_onesimulation_modelchoice(
    onevariation = ls_3markers,
    model = "testDS_lmm")
  # Test output
  # Is it a list?
  expect_type(res_lmm, "list")
  # Does it have the right structure?
  expect_identical(names(res_lmm),
                   c("model", "marker_id", "p_val", "p_adj", "truth"))
})

# Test function_apply_modelcomputations_modelchoice() ---------------------------------------------------------------------------

# Test
test_that("Wrapper to run models", {
  # Same output as in function_to_compute_model_computation_onesimulation_modelchoice()
  # as it is a wrapper around this latest function
  # Set seed
  set.seed(123)
  # Run model - cytoglmm
  res_cytoglmm <- function_apply_modelcomputations_modelchoice(list_combined_output = ls_3markers,
                                                               model = "cytoglmm")

  # Test output
  # Is it a list?
  expect_type(res_cytoglmm, "list")
  # Does it have the right structure?
  expect_identical(names(res_cytoglmm),
                   c("model", "marker_id", "p_val", "p_adj", "truth"))

  # Run model - testDS_limma_random
  res_limma_random <- function_apply_modelcomputations_modelchoice(list_combined_output = ls_3markers,
                                                                   model = "testDS_limma_random")
  # Test output
  # Is it a list?
  expect_type(res_limma_random, "list")
  # Does it have the right structure?
  expect_identical(names(res_limma_random),
                   c("model", "marker_id", "p_val", "p_adj", "truth"))

  # Run model - testDS_limma_fixed
  res_limma_fixed <- function_apply_modelcomputations_modelchoice(list_combined_output = ls_3markers,
                                                                  model = "testDS_limma_fixed")
  # Test output
  # Is it a list?
  expect_type(res_limma_fixed, "list")
  # Does it have the right structure?
  expect_identical(names(res_limma_fixed),
                   c("model", "marker_id", "p_val", "p_adj", "truth"))

  # Run model - testDS_lmm
  res_lmm <- function_apply_modelcomputations_modelchoice(list_combined_output = ls_3markers,
                                                          model = "testDS_lmm")
  # Test output
  # Is it a list?
  expect_type(res_lmm, "list")
  # Does it have the right structure?
  expect_identical(names(res_lmm),
                   c("model", "marker_id", "p_val", "p_adj", "truth"))
})
