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

# Test compute_pwr() ---------------------------------------------------------------------------

# Run model - cytoglmm
res_cytoglmm <- lapply(seq_len(3), function(sim_id){
  res_mod <- function_apply_modelcomputations_modelchoice(list_combined_output = ls_3markers,
                                                          model = "cytoglmm")
})
# Concatenate
df_res_models <- bind_rows(res_cytoglmm, .id = "i")

# Test
test_that("Compute the power", {
  # Set seed
  set.seed(123)
  # Compute power
  pwr <- compute_pwr(df_res_models)

  # Test output
  # Is it a list?
  expect_type(pwr, "list")
  # Does it have the right structure?
  expect_identical(names(pwr),
                   c("model", "marker_id", "power"))
})

# Test compute_variance() ---------------------------------------------------------------------------

# Variance
raw_data_lg <- lapply(seq_len(3), function(sim_id){
  # Generate a single dataset
  sim_data <- function_apply_onesimulation_withmarkerinfo(variation)
  # Observed variance/effect size
  # Merge counts data
  raw_data <- cbind(sim_data$ls_mock_raw_data, "i" = sim_id)
  # Long format
  tidyr::pivot_longer(raw_data,
                      cols = -c("i", "group_id", "donor_id", "sample_id"),
                      names_to = "markers",
                      values_to = "raw_intensity")
})
# Concatenate
df_obs <- bind_rows(raw_data_lg, .id = "i")

# Test
test_that("Compute variance", {
  # Set seed
  set.seed(123)
  # Variance
  obs_variance <- compute_variance(df_obs)

  # Test output
  # Is it a list?
  expect_type(obs_variance, "list")
  # Does it have the right structure?
  expect_identical(names(obs_variance),
                   c("donor_id", "variance"))
})

# Test compute_effectsize() ---------------------------------------------------------------------------

# Test
test_that("Compute variance", {
  # Set seed
  set.seed(123)
  # Effect size
  obs_effectsize <- compute_effectsize(df_obs)

  # Test output
  # Is it a list?
  expect_type(obs_effectsize, "list")
  # Does it have the right structure?
  expect_identical(names(obs_effectsize),
                   c("markers", "effect_size", "observed_FC"))
})
