### Function to run models from the CytoGLMM package ###

#' Run cytoglmm
#'
#' @details Function to run GLMM model from CytoGLMM package.
#' @param mock_dataset data.frame, cell values for each marker.
#'
#' @return list of 3 slots:
#'     - model_fit: fit of the model;
#'     - plot: plot of the effects;
#'     - results_summary: data.frame containing the results of the models for
#'     each marker.
function_run_cytoGLMM <- function(mock_dataset){
  # # Data must be paired
  # if(function_is_data_paired(mock_dataset) == FALSE){
  #   stop("The data is not paired, GLMM cannot be used (please used the bootstrap version)")
  # }
  # Get marker names
  markers_names <- function_extract_marker_names(mock_dataset)
  #Run GLMM
  glmm_fit <- CytoGLMM::cytoglmm(mock_dataset,
                                 protein_names = markers_names,
                                 condition = "group_id", group = "donor_id")
  # Plot the effects
  plot_effects <- plot(glmm_fit)
  # Summary
  res_glmm_fit <- summary(glmm_fit)
  # List of results
  ls_res <- list("model_fit" = glmm_fit,
                 "plot" = plot_effects,
                 "result_summary" = res_glmm_fit)

  return(ls_res)
}

#' Run cytoglm.
#'
#' @details Function to run the Generalized Linear Model with Bootstrap,
#' from the CytoGLMM package.
#' @param mock_dataset data.frame, cell values for each marker.
#' @param nb_bootstrap numeric, number of bootstrap (by defaults nb_bootstrap = 1000).
#'
#' @return list of 3 slots:
#'     - model_fit: fit of the model;
#'     - plot: plot of the effects;
#'     - results_summary: data.frame containing the results of the models for
#'     each marker.
function_run_bootstrapcytoGLMM <- function(mock_dataset, nb_bootstrap = 500){
  # # Data must be paired
  # if(function_is_data_paired(mock_dataset$df_info)){
  #   message("The data is paired, GLMM can also be used")
  # }
  # Get marker names
  markers_names <- function_extract_marker_names(mock_dataset)
  # Run GLM
  glm_fit_bootstrap <- CytoGLMM::cytoglm(mock_dataset,
                                         num_boot = nb_bootstrap,
                                         protein_names = markers_names,
                                         condition = "group_id",
                                         group = "donor_id")
  # Plot the effects
  plot_effects_bootstrap <- plot(glm_fit_bootstrap)
  # Summary
  res_glm_fit_bootstrap <- summary(glm_fit_bootstrap)

  # List of results
  ls_res <- list("model_fit" = glm_fit_bootstrap,
                 "plot" = plot_effects_bootstrap,
                 "result_summary" = res_glm_fit_bootstrap)

  return(ls_res)
}

### Functions to run the models from the diffcyt package ###

#' Compute cells counts and medians.
#'
#' @details Function to calculate features for diffcyt package models.
#'
#' @param mock_flowset data.frame, cell values for each marker.
#'
#' @return list with cell counts and medians for each markers.
function_compute_diffcyt_features <- function(mock_flowset){
  # Calculate features
  # Cluster cell counts
  d_mock_flowset_counts <- diffcyt::calcCounts(mock_flowset)
  # Cluster medians
  d_mock_flowset_medians <- diffcyt::calcMedians(mock_flowset)
  # Return
  return(list("counts" = d_mock_flowset_counts,
              "medians" = d_mock_flowset_medians))
}

#' Design and contrast matrices for diffcyt-DS-limma with random effect.
#'
#' @details Function to create the design matrix and contrast for the
#' diffcyt-DS-limma model with random effect from the diffcyt package.
#'
#' @param df_experiment_info data.frame, information about the experiment (donor IDs,
#' group IDs and sample IDs).
#'
#' @return list with 3 slots:
#'     - design_matrix: the design matrix;
#'     - contrast: the contrast matrix;
#'     - effect: the specification of the "random" effect.
function_desigmat_contrast_diffcytDSlimma_randomeffect <- function(df_experiment_info){
  # Create design matrix
  designmatrix <- diffcyt::createDesignMatrix(df_experiment_info,
                                              cols_design = c("group_id"))
  # Create contrast
  contrast <- diffcyt::createContrast(c(0, 1))
  # Assign a variable for the effect
  effect <- "random"
  # Return
  return(list("design_matrix" = designmatrix,
              "contrast" = contrast,
              "effect" = effect))
}

#' Design and contrast matrices for diffcyt-DS-limma with fixed effect.
#'
#' @details Function to create the design matrix and contrast for the diffcyt-DS-limma model
#'with fixed effect.
#' @param df_experiment_info data.frame, information about the experiment (donor IDs,
#' group IDs and sample IDs).
#'
#' @return list with 3 slots:
#'     - design_matrix: the design matrix;
#'     - contrast: the contrast matrix;
#'     - effect: the specification of the "fixed" effect.
function_desigmat_contrast_diffcytDSlimma_fixedeffect <- function(df_experiment_info){
  # Create design matrix
  designmatrix <- diffcyt::createDesignMatrix(df_experiment_info,
                                              cols_design = c("group_id", "donor_id"))
  # Create contrast
  contrast <- diffcyt::createContrast(c(0, 1, rep(0, dim(designmatrix)[2]-2)))
  # Assign a variable for the effect
  effect <- "fixed"
  # Return
  return(list("design_matrix" = designmatrix,
              "contrast" = contrast,
              "effect" = effect))
}

#' Run diffcyt-DS-limma model.
#'
#' @details Function to run diffcyt-DS-limma model from the diffcyt package.
#'
#' @param ls_desigmat_contrast list, with design and contrast matrices.
#' @param df_experiment_info data.frame, information about the experiment (donor IDs,
#' group IDs and sample IDs).
#' @param ls_features list, with cell counts and medians for each markers.
#'
#' @return list with 2 slots:
#'     - model_fit: the model fit;
#'     - result_summary: results of the model for each marker.
function_run_diffcytDSlimma <- function(ls_desigmat_contrast,
                                        df_experiment_info,
                                        ls_features){
  # Run the model depending on the effects (random or fixed)
  if(ls_desigmat_contrast$effect == "random"){
    message("Running diffcyt-DS-limma model with random effect")
    res_mock_flowset <- diffcyt::testDS_limma(
      d_counts = ls_features$counts,
      d_medians = ls_features$medians,
      design = ls_desigmat_contrast$design_matrix,
      contrast = ls_desigmat_contrast$contrast,
      block_id = df_experiment_info$donor_id)
  } else if(ls_desigmat_contrast$effect == "fixed"){
    message("Running diffcyt-DS-limma model with fixed effect")
    res_mock_flowset <- diffcyt::testDS_limma(
      d_counts = ls_features$counts,
      d_medians = ls_features$medians,
      design = ls_desigmat_contrast$design_matrix,
      contrast = ls_desigmat_contrast$contrast)
  }
  # Results
  table_res_mock_flowset <- SummarizedExperiment::rowData(res_mock_flowset, format_vals = TRUE)
  # Return
  ls_res <- list("model_fit" = res_mock_flowset, "result_summary" = table_res_mock_flowset)
  return(ls_res)
}

#' Formula and contrast matrix for diffcyt-DS-LMM with random effect.
#'
#' @details Function to create formula and contrast for diffcyt-DS-LMM with random effect.
#'
#' @param df_experiment_info data.frame, information about the experiment (donor IDs,
#' group IDs and sample IDs).
#'
#' @return list with 3 slots:
#'     - formula: the formula for the model;
#'     - contrast: the contrast matrix.
function_formula_contrast_diffcytDSLMM_randomeffect <- function(df_experiment_info){
  # Create formula
  formula_withrandomeffect <-
    diffcyt::createFormula(df_experiment_info,
                           cols_fixed = c("group_id"),
                           cols_random = c("donor_id"))
  # Create contrast
  contrast_withrandomeffect <- diffcyt::createContrast(c(0, 1)) # no 0 for the random effect
  # Return
  return(list("formula" = formula_withrandomeffect,
              "contrast" = contrast_withrandomeffect))
}

#' Run diffcyt-DS-LMM model.
#'
#' @details Function to run diffcyt-DS-LMM model from the diffcyt package.
#'
#' @param ls_form_contrast list, with formula and constrat matrix.
#' @param df_experiment_info data.frame, information about the experiment (donor IDs,
#' group IDs and sample IDs).
#' @param ls_features list, with cell counts and medians for each markers.
#'
#' @return list with 2 slots:
#'     - model_fit: the model fit;
#'     - result_summary: results of the model for each marker.
function_run_diffcytDSLMM <- function(ls_form_contrast, df_experiment_info, ls_features){
  # DS test within clusters with LMM
  res_mock_flowset <- diffcyt::testDS_LMM(
    d_counts = ls_features$counts,
    d_medians = ls_features$medians,
    formula = ls_form_contrast$formula,
    contrast = ls_form_contrast$contrast)
  # Results
  table_res_mock_flowset <- SummarizedExperiment::rowData(res_mock_flowset,
                                                          format_vals = TRUE)
  # Return
  ls_res <- list("model_fit" = res_mock_flowset,
                 "result_summary" = table_res_mock_flowset)
  return(ls_res)
}

#' Run diffcyt pipeline
#'
#' @details Function to run diffcyt pipeline. We do not used the diffcyt function
#' directly because of the limma model run differently with the effects.
#'
#' @param onevariation list, of simulated data (output of the
#'   function_wrapper_apply_simulation_nbtimes function).
#' @param model character, model to run: "limma" or "LMM".
#' @param effect character, effect: "random" or "fixed".
#'
#' @return list with 2 slots:
#'     - model_fit: the model fit;
#'     - result_summary: results of the model for each marker.
function_run_diffcyt_full_pipeline <- function(onevariation,
                                               model = c("limma", "LMM"),
                                               effect = c("random", "fixed")){
  # Prepare data
  # Mock dataset
  mock_dataset <- onevariation$ls_mock_data
  # Experimental info
  df_experiment_info <- onevariation$df_info
  # Marker names
  markers_names <- function_extract_marker_names(mock_dataset)
  # Create metadata
  # colData
  mock_flowset_markerinfo <- data.frame("channel_name" = markers_names,
                                        "marker_name" = markers_names,
                                        "marker_class" = "state")
  # rowData
  se_df_info <- data.frame(mock_dataset[,c("group_id", "donor_id", "sample_id")],
                           "cluster_id" = "pop1",
                           stringsAsFactors = TRUE)
  # Number of cells per sample
  se_n_cells <- rep(unique(onevariation$variation$nb_cell_per_sample),
                    each = length(onevariation$df_info$sample_id))
  names(se_n_cells) <- paste0("Sample",
                              seq_len(length(onevariation$df_info$sample_id)))
  # metadata
  met <- list("experiment_info" = data.frame(onevariation$df_info,
                                             "cluster_id" = "pop1"),
              "n_cells" = se_n_cells)
  # Create SummarizedExperiment
  d_flowset <- SummarizedExperiment::SummarizedExperiment(
    assays = list("exprs" = as.matrix(mock_dataset[, markers_names])),
    rowData = se_df_info,
    colData = mock_flowset_markerinfo,
    metadata = met)
  # Compute features
  ls_features <- function_compute_diffcyt_features(d_flowset)
  # Create design matrix or formula and constrast depending on choosen model
  if(model == "LMM"){
    message("Run the LMM model, with random effect")
    # Create formula and contrast
    ls_form <- function_formula_contrast_diffcytDSLMM_randomeffect(df_experiment_info)
    # Run the model
    res <- function_run_diffcytDSLMM(ls_form, df_experiment_info, ls_features)
  } else if(model == "limma" & effect == "random"){
    message("Run the limma model with random effect")
    # Create design matrix and contrast
    ls_desig <- function_desigmat_contrast_diffcytDSlimma_randomeffect(df_experiment_info)
    # Run the model
    res <- function_run_diffcytDSlimma(ls_desig, df_experiment_info, ls_features)
  } else if(model == "limma" & effect == "fixed"){
    message("Run the limma model with fixed effect")
    # Create design matrix and contrast
    ls_desig <- function_desigmat_contrast_diffcytDSlimma_fixedeffect(df_experiment_info)
    # Run the model
    res <- function_run_diffcytDSlimma(ls_desig, df_experiment_info, ls_features)
  } else {
    stop("Model not correctly specified")
  }
  # Return
  return(res)
}





