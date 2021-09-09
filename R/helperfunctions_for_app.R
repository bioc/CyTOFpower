#' Compute power
#'
#' @details Compute the power based on the model values
#'
#' @param model_values data.frame, output of run_models()
#' @param alpha numeric, significance level
compute_pwr <- function(model_values,
                        alpha = 0.05){
  model_values %>%
    dplyr::filter(.data$truth == 1) %>% # Filter the values with a true fold change
    dplyr::group_by(.data$model, .data$marker_id) %>% # Group by models and marker ID
    dplyr::summarize(power = mean(.data$p_adj < alpha)) # Compute the power
}

#' Compute variance
#'
#' @details Compute the observed variance in the data
#'
#' @param raw_data_lg data.frame, cells values in long format
compute_variance <- function(raw_data_lg){
  # Sample variance
  obs_sample_variance <- raw_data_lg %>%
    dplyr::group_by(.data$i, .data$donor_id) %>%
    dplyr::summarize(variance = stats::var(.data$raw_intensity))
  # Mean of the sample variance
  mean_obs_sample_variance <- obs_sample_variance %>%
    dplyr::group_by(.data$donor_id) %>%
    dplyr::summarize(variance = round(mean(.data$variance), digits = 2))

  return(mean_obs_sample_variance)
}

#' Compute effect size
#'
#' @details Compute observed Cohen's effect size and observed fold change
#'
#' @param raw_data_lg data.frame, cells values in long format
compute_effectsize <- function(raw_data_lg){
  # Effect size
  ## Mean
  mean_raw_data_summary_1 <- raw_data_lg %>%
    dplyr::group_by(.data$group_id, .data$markers) %>%
    dplyr::summarize(mean = mean(.data$raw_intensity))
  ## Mean difference
  mean_raw_data_summary_1_wd <- mean_raw_data_summary_1 %>%
    tidyr::pivot_wider(id_cols = .data$markers,
                       names_from = .data$group_id,
                       values_from = .data$mean)
  mean_diff_raw_data_summary_1_wd <- dplyr::summarize(mean_raw_data_summary_1_wd,
                                                      mean_diff = .data$B - .data$A,
                                                      .data$markers)
  ## Standard deviation
  sd_raw_data_summary_1 <- raw_data_lg %>%
    dplyr::group_by(.data$i, .data$markers)  %>%
    dplyr::summarize(standard_deviation = stats::sd(.data$raw_intensity))
  sd_mean_raw_data_summary_1 <- sd_raw_data_summary_1 %>%
    dplyr::group_by(.data$markers)  %>%
    dplyr::summarize(mean_sd = mean(.data$standard_deviation))
  ## Combine
  cb_raw_data <- dplyr::left_join(mean_diff_raw_data_summary_1_wd, sd_mean_raw_data_summary_1)
  eff_size_raw_data <- dplyr::summarize(cb_raw_data,
                                        effect_size = round(.data$mean_diff / .data$mean_sd, digits = 1),
                                        .data$markers)

  # Fold Change
  obs_fc <- dplyr::summarize(mean_raw_data_summary_1_wd,
                             observed_FC = round(.data$B / .data$A, digits = 1),
                             .data$markers)

  # Return
  res <- dplyr::left_join(eff_size_raw_data, obs_fc)
  return(res[,c("markers", "effect_size", "observed_FC")])
}
