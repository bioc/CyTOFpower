# Functions to summarize results from the different models

#' Summarize data and results.
#'
#' @details Function to do a summary of the tested data and model's results
#'   for the CytoGLMM and diffcyt packages.
#'
#' @param summary_from_model list, output from the functions running the models.
#' @param package character, package used to run the test: "CytoGLMM" or "diffcyt".
#'
#' @return data.frame of results for each simulation.
function_summary_results_models <- function(summary_from_model, package){

  # Get package argument
  package <- match.arg(package, choices = c("CytoGLMM", "diffcyt"))
  message("Using package: ", package)

  ## Rename colnames
  ### if the cytoGLMM package has been used, the colnames in the results table is "protein_name"
  if(package == "CytoGLMM"){
    res_sum <- as.data.frame(summary_from_model$result_summary)
    colnames(res_sum) <- c("marker_id", "p_val", "p_adj")
  } else if(package =="diffcyt"){ ## if the diffcyt package has been used, the colnames in the results table is "marker_id"
    res_sum <- as.data.frame(summary_from_model$result_summary[,c("marker_id", "p_val", "p_adj")])
  }
  # Return
  return(as.data.frame(res_sum))
}

#' Run DS tests.
#'
#' @details Function to run the DS tests through the different models.
#'
#' @param simulation numeric, counter for each simulation.
#' @param onevariation list, of simulated data (output of the
#'   function_wrapper_apply_simulation_nbtimes function).
#'
#' @return data.frame of results for each simulation all models combined.
function_to_compute_model_computation_from_simulations <- function(simulation, onevariation){
  # Split info from the mock datasets
  # nb_donor <- onevariation$variation$nb_donor
  # delta <- onevariation$variation$delta
  # subject_effect <- onevariation$variation$subject_effect
  # nb_DE_marker <- onevariation$variation$nb_DE_marker
  df_exp_info <- onevariation$df_info
  onemockdataset <- onevariation$ls_mock_data[[simulation]]

  # Names DE markers
  #vec_names_DEmarkers <- function_extract_DE_marker_names(data = onemockdataset, nb_DE_marker = nb_DE_marker)
  #names_DEmarkers <- paste(vec_names_DEmarkers, collapse = ";")
  vec_names_DEmarkers <- onevariation$DE_markers_names
  names_DEmarkers <- paste(vec_names_DEmarkers, collapse = ";")

  # CytoGLMM
  ## GLMM
  ls_glmm <- function_run_cytoGLMM(onemockdataset)
  #### All p-values
  allres_glmm <- function_summary_results_models(ls_glmm, package = "CytoGLMM")
  # Bootstrap
  ls_bootstrap_glm <- function_run_bootstrapcytoGLMM(onemockdataset, nb_bootstrap = 50)
  allres_glm <- function_summary_results_models(ls_bootstrap_glm, package = "CytoGLMM")

  # Diffcyt
  ## Run models
  ### limma random
  ls_diffcyt_limma_random <- function_run_diffcyt_full_pipeline(df_experiment_info = df_exp_info,
                                                                mock_dataset = onemockdataset,
                                                                model = "limma",
                                                                effect = "random")
  #### All p-values
  allres_diffcyt_limma_random <- function_summary_results_models(ls_diffcyt_limma_random, package = "diffcyt")
  ### limma fixed
  ls_diffcyt_limma_fixed <- function_run_diffcyt_full_pipeline(df_experiment_info = df_exp_info,
                                                               mock_dataset = onemockdataset,
                                                               model = "limma",
                                                               effect = "fixed")
  #### All p-values
  allres_diffcyt_limma_fixed <- function_summary_results_models(ls_diffcyt_limma_fixed, package = "diffcyt")
  ### LMM
  ls_diffcyt_LMM_random <- function_run_diffcyt_full_pipeline(df_experiment_info = df_exp_info,
                                                              mock_dataset = onemockdataset,
                                                              model = "LMM",
                                                              effect = "random")
  #### All p-values
  allres_diffcyt_LMM_random <- function_summary_results_models(ls_diffcyt_LMM_random, package = "diffcyt")

  # Combine p-val tables
  allpval_res <- dplyr::bind_rows("glmm" = allres_glmm,
                                  "glm" = allres_glm,
                                  "limma_random" = allres_diffcyt_limma_random,
                                  "limma_fixed" = allres_diffcyt_limma_fixed,
                                  "lmm_random" = allres_diffcyt_LMM_random,
                                  .id = "model")
  # allpval_res <- cbind(allpval_res, data.frame("nb_donor" = nb_donor,
  #                                              "delta" = delta,
  #                                              "subject_effect" = subject_effect,
  #                                              "nb_DE_marker" = nb_DE_marker,
  #                                              "names_DEmarkers" = names_DEmarkers))
  allpval_res <- cbind(allpval_res, onevariation$variation)
  # Add the information about the truth, i.e. which markers have DE by design
  #allpval_res$truth <- 0
  allpval_res <- dplyr::mutate(allpval_res,
                               truth = ifelse(.data$marker_id %in% as.character(vec_names_DEmarkers), 1, 0))
  #allpval_res %>% dplyr::mutate(truth = replace(truth, marker_id == as.character(vec_names_DEmarkers), 1))

  # Return
  return(allpval_res)
}

#' Run DS tests.
#'
#' @details Function to run the DS tests through the different models.
#'
#' @param simulation numeric, counter for each simulation.
#' @param onevariation list, of simulated data (output of the
#'   function_wrapper_apply_simulation_nbtimes function).
#' @param model vector, name(s) of models to test.
#'
#' @return data.frame of results for each simulation all models combined.
function_to_compute_model_computation_from_simulations_modelchoice <- function(simulation, onevariation, model){
  # Split info from the mock datasets
  df_exp_info <- onevariation$df_info
  onemockdataset <- onevariation$ls_mock_data[[simulation]]

  # Names DE markers
  vec_names_DEmarkers <- onevariation$DE_markers_names
  names_DEmarkers <- paste(vec_names_DEmarkers, collapse = ";")

  # Initialize results list
  ls_res <- list()

  # cytoGLMM
  if (any(model %in% "cytoglmm")) {
    # CytoGLMM - GLMM
    ls_glmm <- function_run_cytoGLMM(onemockdataset)
    #### All p-values
    allres_glmm <- function_summary_results_models(ls_glmm, package = "CytoGLMM")
    ls_res[["cytoglmm"]] <- allres_glmm
  }
  if (any(model %in% "cytoglm")) {
    # CytoGLMM - Bootstrap
    ls_bootstrap_glm <- function_run_bootstrapcytoGLMM(onemockdataset, nb_bootstrap = 50)
    allres_glm <- function_summary_results_models(ls_bootstrap_glm, package = "CytoGLMM")
    ls_res[["cytoglm"]] <- allres_glm
  }

  # Diffcyt
  ## Run models
  if (any(model %in% "testDS_limma_random")) {
    ### limma random
    ls_diffcyt_limma_random <- function_run_diffcyt_full_pipeline(df_experiment_info = df_exp_info,
                                                                  mock_dataset = onemockdataset,
                                                                  model = "limma",
                                                                  effect = "random")
    #### All p-values
    allres_diffcyt_limma_random <- function_summary_results_models(ls_diffcyt_limma_random, package = "diffcyt")
    ls_res[["testDS_limma_random"]] <- allres_diffcyt_limma_random
  }

  if (any(model %in% "testDS_limma_fixed")) {
    ### limma fixed
    ls_diffcyt_limma_fixed <- function_run_diffcyt_full_pipeline(df_experiment_info = df_exp_info,
                                                                 mock_dataset = onemockdataset,
                                                                 model = "limma",
                                                                 effect = "fixed")
    #### All p-values
    allres_diffcyt_limma_fixed <- function_summary_results_models(ls_diffcyt_limma_fixed, package = "diffcyt")
    ls_res[["testDS_limma_fixed"]] <- allres_diffcyt_limma_fixed
  }

  if (any(model %in% "testDS_lmm") ) {
    ### LMM
    ls_diffcyt_LMM_random <- function_run_diffcyt_full_pipeline(df_experiment_info = df_exp_info,
                                                                mock_dataset = onemockdataset,
                                                                model = "LMM",
                                                                effect = "random")
    #### All p-values
    allres_diffcyt_LMM_random <- function_summary_results_models(ls_diffcyt_LMM_random, package = "diffcyt")
    ls_res[["testDS_lmm"]] <- allres_diffcyt_LMM_random
  }

  # Combine p-val tables
  allpval_res <- bind_rows(ls_res, .id = "model")
  # allpval_res <- cbind(allpval_res, onevariation$variation)
  # Add the information about the truth, i.e. which markers have DE by design
  allpval_res <- dplyr::mutate(allpval_res,
                               truth = ifelse(.data$marker_id %in% as.character(vec_names_DEmarkers), 1, 0))

  # Return
  return(allpval_res)
}
