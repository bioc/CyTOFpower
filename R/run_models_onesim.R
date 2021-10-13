#' Run DS tests for one simulation.
#'
#' @details Wrapper to run the models through the different simulations.
#'
#' @param list_combined_output list, of simulated data (output of the
#' function_apply_onesimulation_withmarkerinfo function).
#' @param model vector, name(s) of models to test.
#'
#' @return data.frame of results for each simulation.
function_apply_modelcomputations_modelchoice <- function(list_combined_output,
                                                          model = c("cytoglmm",
                                                                    "cytoglm",
                                                                    "testDS_limma_random",
                                                                    "testDS_limma_fixed",
                                                                    "testDS_lmm")) {

  print(list_combined_output$variation)
  # Repeat model computation
  res_output <-
    function_to_compute_model_computation_onesimulation_modelchoice(onevariation = list_combined_output,
                                                                    model = model)
  row.names(res_output) <- seq_len(dim(res_output)[1])
  # Return
  return(res_output)
}

#' Run DS tests for one simulation.
#'
#' @details Function to run the DS tests through the different models.
#'
#' @param onevariation list, of simulated data (output of the
#'   function_wrapper_apply_simulation_nbtimes function).
#' @param model vector, name(s) of models to test.
#'
#' @return data.frame of results for each simulation all models combined.
function_to_compute_model_computation_onesimulation_modelchoice <- function(onevariation,
                                                                            model){
  # Split info from the mock datasets
  onemockdataset <- onevariation$ls_mock_data

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
    ls_diffcyt_limma_random <- function_run_diffcyt_full_pipeline(onevariation = onevariation,
                                                                  model = "limma",
                                                                  effect = "random")
    #### All p-values
    allres_diffcyt_limma_random <- function_summary_results_models(ls_diffcyt_limma_random, package = "diffcyt")
    ls_res[["testDS_limma_random"]] <- allres_diffcyt_limma_random
  }

  if (any(model %in% "testDS_limma_fixed")) {
    ### limma fixed
    ls_diffcyt_limma_fixed <- function_run_diffcyt_full_pipeline(onevariation = onevariation,
                                                                 model = "limma",
                                                                 effect = "fixed")
    #### All p-values
    allres_diffcyt_limma_fixed <- function_summary_results_models(ls_diffcyt_limma_fixed, package = "diffcyt")
    ls_res[["testDS_limma_fixed"]] <- allres_diffcyt_limma_fixed
  }

  if (any(model %in% "testDS_lmm") ) {
    ### LMM
    ls_diffcyt_LMM_random <- function_run_diffcyt_full_pipeline(onevariation = onevariation,
                                                                model = "LMM",
                                                                effect = "random")
    #### All p-values
    allres_diffcyt_LMM_random <- function_summary_results_models(ls_diffcyt_LMM_random, package = "diffcyt")
    ls_res[["testDS_lmm"]] <- allres_diffcyt_LMM_random
  }

  # Combine p-val tables
  allpval_res <- dplyr::bind_rows(ls_res, .id = "model")
  # Add the information about the truth, i.e. which markers have DE by design
  allpval_res <- dplyr::mutate(allpval_res,
                               truth = ifelse(.data$marker_id %in% as.character(vec_names_DEmarkers), 1, 0))

  # Return
  return(allpval_res)
}
