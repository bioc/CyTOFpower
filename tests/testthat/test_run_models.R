# Test function to run the models

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

# Create a summarized experiment with the mock dataset
# Marker names
markers_names <- function_extract_marker_names(ls_3markers$ls_mock_data)
# Create metadata
# colData
mock_flowset_markerinfo <- data.frame("channel_name" = markers_names,
                                      "marker_name" = markers_names,
                                      "marker_class" = "state")
# rowData
se_df_info <- data.frame(ls_3markers$ls_mock_data[,c("group_id", "donor_id", "sample_id")],
                         "cluster_id" = "pop1",
                         stringsAsFactors = TRUE)
# Number of cells per sample
se_n_cells <- rep(unique(ls_3markers$variation$nb_cell_per_sample),
                  each = length(ls_3markers$df_info$sample_id))
names(se_n_cells) <- paste0("Sample",
                            seq_len(length(ls_3markers$df_info$sample_id)))
# metadata
met <- list("experiment_info" = data.frame(ls_3markers$df_info,
                                           "cluster_id" = "pop1"),
            "n_cells" = se_n_cells)
# Create SummarizedExperiment
d_sumexp <- SummarizedExperiment::SummarizedExperiment(
  assays = list("exprs" = as.matrix(ls_3markers$ls_mock_data[, markers_names])),
  rowData = se_df_info,
  colData = mock_flowset_markerinfo,
  metadata = met)


# test function_run_cytoGLMM() ---------------------------------------------------------------------------

# Expected p-values data.frame
df_res_sum <- tibble::tibble(
  "protein_name" = c("m2", "m3", "m1"),
  "pvalues_unadj" = c("m2" = 0.408, "m3" = 0.760, "m1" = 0.866),
  "pvalues_adj" = c("m2" = 0.866, "m3" = 0.866, "m1" = 0.866))

# Test
test_that("Run the CytoGLMM - GLMM model", {
  # Set seed
  set.seed(123)
  # Run model
  expect_warning(cytoglmm_res <- function_run_cytoGLMM(
    mock_dataset = ls_3markers$ls_mock_data))
  # Test output
  # Is it a list?
  expect_type(cytoglmm_res, "list")
  # Does it have the right structure?
  expect_identical(names(cytoglmm_res),
                   c("model_fit", "plot", "result_summary"))
  # Are the p-values correct?
  expect_equal(cytoglmm_res$result_summary,
               df_res_sum,
               tolerance = 1e-3)
})

# test function_compute_diffcyt_features() ---------------------------------------------------------------------------

# Expected counts in the assays
v_exp_counts <- matrix(rep(3, 6),
                       nrow = 1,
                       dimnames = list(c("pop1"),
                                       c("sample_1",
                                         "sample_4",
                                         "sample_2",
                                         "sample_5",
                                         "sample_3",
                                         "sample_6")))
# Expected medians in the assays
v_exp_medians <- matrix(c(1.529, 0.732, 1.350, 1.350, 1.753, 1.350),
                        nrow = 1,
                        dimnames = list(c("pop1"),
                                        c("sample_1",
                                          "sample_4",
                                          "sample_2",
                                          "sample_5",
                                          "sample_3",
                                          "sample_6")))


# Test
test_that("Compute the diffcyt features", {
  # Set seed
  set.seed(123)
  # Compute the features
  ls_features <- function_compute_diffcyt_features(d_sumexp)
  # Test output
  # Is it a list?
  expect_type(ls_features, "list")
  # Does it have the right structure?
  expect_identical(names(ls_features),
                   c("counts", "medians"))
  # Are the values correct?
  # Counts
  expect_equal(SummarizedExperiment::assay(ls_features$counts),
               v_exp_counts)
  # Median
  expect_equal(SummarizedExperiment::assay(ls_features$medians),
               v_exp_medians,
               tolerance = 1e-3)
})

# test function_desigmat_contrast_diffcytDSlimma_randomeffect() ---------------------------------------------------------------------------

# Expected design matrix
mat_design_mat_randomeffect <- matrix(c(rep(1, 6), rep(c(0, 1), times = 3)),
                                      ncol = 2,
                                      dimnames = list(1:6,
                                                      c("(Intercept)", "group_idB")))

# Test
test_that("Contrast for limma model with random effect", {
  # Set seed
  set.seed(123)
  # Generate the contrast
  ls_contrast <- function_desigmat_contrast_diffcytDSlimma_randomeffect(df_exp_info)
  # Test output
  # Is it a list?
  expect_type(ls_contrast, "list")
  # Does it have the right structure?
  expect_identical(names(ls_contrast),
                   c("design_matrix", "contrast", "effect"))
  # Are the values of the matrix equal to the expected values?
  expect_equal(ls_contrast$design_matrix[,"(Intercept)"],
               mat_design_mat_randomeffect[,"(Intercept)"])
  expect_equal(ls_contrast$design_matrix[,"group_idB"],
               mat_design_mat_randomeffect[,"group_idB"])
  expect_identical(ls_contrast$contrast,
                   matrix(c(0, 1), ncol = 1))
  expect_identical(ls_contrast$effect,
                   "random")
})

# test function_desigmat_contrast_diffcytDSlimma_fixedeffect() ---------------------------------------------------------------------------

# Expected design matrix
mat_design_mat_donorids <-
  matrix(c(0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1),
         ncol = 2,
         dimnames = list(1:6,
                         c("donor_iddonor_2", "donor_iddonor_3")))
mat_design_mat_fixedeffect <- cbind(mat_design_mat_randomeffect,
                                    mat_design_mat_donorids)

# Test
test_that("Contrast for limma model with fixed effect", {
  # Set seed
  set.seed(123)
  # Generate the contrast
  ls_contrast <- function_desigmat_contrast_diffcytDSlimma_fixedeffect(df_exp_info)
  # Test output
  # Is it a list?
  expect_type(ls_contrast, "list")
  # Does it have the right structure?
  expect_identical(names(ls_contrast),
                   c("design_matrix", "contrast", "effect"))
  # Are the values of the matrix equal to the expected values?
  expect_equal(ls_contrast$design_matrix[,"(Intercept)"],
               mat_design_mat_fixedeffect[,"(Intercept)"])
  expect_equal(ls_contrast$design_matrix[,"group_idB"],
               mat_design_mat_fixedeffect[,"group_idB"])
  expect_equal(ls_contrast$design_matrix[,"donor_iddonor_2"],
               mat_design_mat_fixedeffect[,"donor_iddonor_2"])
  expect_equal(ls_contrast$design_matrix[,"donor_iddonor_3"],
               mat_design_mat_fixedeffect[,"donor_iddonor_3"])
  expect_identical(ls_contrast$contrast,
                   matrix(c(0, 1, 0, 0), ncol = 1))
  expect_identical(ls_contrast$effect,
                   "fixed")
})

# test function_run_diffcytDSlimma() ---------------------------------------------------------------------------

# Test
test_that("Run diffcyt-DS-limma model", {
  # Set seed
  set.seed(123)
  # Run the model
  # Compute the features
  ls_features <- function_compute_diffcyt_features(d_sumexp)

  # With fixed effect
  # Compute the contrast
  ls_desigmat_contrast_f <- function_desigmat_contrast_diffcytDSlimma_fixedeffect(df_exp_info)
  # Run model
  ls_res_model_f <- function_run_diffcytDSlimma(ls_desigmat_contrast_f,
                                                df_experiment_info,
                                                ls_features)
  # Test the output
  # Is it a list?
  expect_type(ls_res_model_f, "list")
  # Does it have the right structure?
  expect_identical(names(ls_res_model_f),
                   c("model_fit", "result_summary"))
  # Are the p-values equal to the expected ones?
  expect_equal(ls_res_model_f$result_summary$p_adj,
               c(0.127, 0.049, 0.737),
               tolerance = 1e-3)

  # With random effect
  # Compute the contrast
  ls_desigmat_contrast_r <- function_desigmat_contrast_diffcytDSlimma_fixedeffect(df_exp_info)
  # Run the model
  ls_res_model_r <- function_run_diffcytDSlimma(ls_desigmat_contrast_r,
                                                df_experiment_info,
                                                ls_features)
  # Test the output
  # Is it a list?
  expect_type(ls_res_model_r, "list")
  # Does it have the right structure?
  expect_identical(names(ls_res_model_r),
                   c("model_fit", "result_summary"))
  # Are the p-values equal to the expected ones?
  expect_equal(ls_res_model_r$result_summary$p_adj,
               c(0.127, 0.049, 0.737),
               tolerance = 1e-3)
})

# test function_formula_contrast_diffcytDSLMM_randomeffect() ---------------------------------------------------------------------------

# Test
test_that("Contrast for LMM model", {
  # Set seed
  set.seed(123)
  # Generate formula and contrast
  ls_contrast_LMM <- function_formula_contrast_diffcytDSLMM_randomeffect(df_exp_info)
  # Test output
  # Is it a list?
  expect_type(ls_contrast_LMM, "list")
  # Does it have the right structure?
  expect_identical(names(ls_contrast_LMM),
                   c("formula", "contrast"))
  # Check the formula
  expect_equal(ls_contrast_LMM$formula$formula,
               formula("y ~ group_id + (1 | donor_id)"),
               ignore_formula_env	= TRUE)
  # Is the random_terms equal TRUE?
  expect_true(ls_contrast_LMM$formula$random_terms)
  # Are the data equal to the experimental information?
  expect_identical(ls_contrast_LMM$formula$data,
                   df_exp_info[, c("group_id", "donor_id")])
  # Check contrast
  expect_identical(ls_contrast_LMM$contrast,
                   matrix(c(0, 1), ncol = 1))
})

# test function_run_diffcytDSLMM() ---------------------------------------------------------------------------

# Test
test_that("Run diffcyt-DS-LMM model", {
  # Set seed
  set.seed(123)
  # Generate formula/contrast
  ls_contrast_LMM <- function_formula_contrast_diffcytDSLMM_randomeffect(
    df_experiment_info = df_exp_info)
  # Compute the features
  ls_features <- function_compute_diffcyt_features(d_sumexp)
  # Run the model
  ls_res_model_lmm <- function_run_diffcytDSLMM(ls_form_contrast = ls_contrast_LMM,
                                                df_experiment_info = df_exp_info,
                                                ls_features = ls_features)

  # Test the output
  # Is it a list?
  expect_type(ls_res_model_lmm, "list")
  # Does it have the right structure?
  expect_identical(names(ls_res_model_lmm),
                   c("model_fit", "result_summary"))
  # Are the p-values equal to the expected ones?
  expect_equal(ls_res_model_lmm$result_summary$p_adj,
               c(0.123, 0, 0.544),
               tolerance = 1e-3)
})

# test function_run_diffcyt_full_pipeline() ---------------------------------------------------------------------------

# Test
test_that("Run the full diffcyt pipeline", {
  # Set seed
  set.seed(123)
  # In this test we check that the switch between the different models is
  # done correctly

  # diffcyt-DS-limma random
  # Check message
  expect_message(function_run_diffcyt_full_pipeline(onevariation = ls_3markers,
                                                    model = c("limma"),
                                                    effect = c("random")),
                 regexp = "Run the limma model with random effect")
  # diffcyt-DS-limma fixed
  # Check message
  expect_message(function_run_diffcyt_full_pipeline(onevariation = ls_3markers,
                                                    model = c("limma"),
                                                    effect = c("fixed")),
                 regexp = "Run the limma model with fixed effect")
  # diffcyt-DS-LMM
  # Check message
  expect_message(function_run_diffcyt_full_pipeline(onevariation = ls_3markers,
                                                    model = c("LMM")),
                 regexp = "Run the LMM model with random effect")
})
