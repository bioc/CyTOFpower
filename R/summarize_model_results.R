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

  # Rename colnames
  # if the cytoGLMM package has been used, the column name in the results table
  # is "protein_name"
  if(package == "CytoGLMM"){
    res_sum <- as.data.frame(summary_from_model$result_summary)
    colnames(res_sum) <- c("marker_id", "p_val", "p_adj")
  } else if(package =="diffcyt"){
    # if the diffcyt package has been used, the column name in the results table
    # is "marker_id"
    res_sum <- as.data.frame(summary_from_model$result_summary[,c("marker_id",
                                                                  "p_val",
                                                                  "p_adj")])
  }
  # Return
  return(as.data.frame(res_sum))
}
