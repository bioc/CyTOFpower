# Compute power
compute_pwr <- function(model_values, # Output of run_models()
                        alpha = 0.05){ # Significance level
  # Compute the power for each marker
  model_values %>%
    dplyr::filter(truth == 1) %>%
    dplyr::group_by(model,marker_id) %>% #  rho, nb_donor, subject_effect, nb_DE_marker,  nb_cell_per_sample, total_nb_marker
    dplyr::summarize(power = mean(p_adj < alpha))
}

# Compute variance
compute_variance <- function(raw_data_lg){
  # Observed sample variance
  obs_sample_variance <- raw_data_lg %>%
    dplyr::group_by(i, donor_id) %>%
    dplyr::summarize(variance = var(raw_intensity))
  mean_obs_sample_variance <- obs_sample_variance %>%
    dplyr::group_by(donor_id) %>%
    dplyr::summarize(mean_variance = round(mean(variance), digits = 2))

  return(mean_obs_sample_variance)
}
# Compute variance
compute_effectsize <- function(raw_data_lg){
  # Effect size
  ## Mean
  mean_raw_data_summary_1 <- raw_data_lg %>%
    dplyr::group_by(group_id, markers) %>%
    dplyr::summarize(mean = mean(raw_intensity))
  ## Mean difference
  mean_raw_data_summary_1_wd <- mean_raw_data_summary_1 %>%
    tidyr::pivot_wider(id_cols = markers, names_from = group_id, values_from = mean)
  mean_diff_raw_data_summary_1_wd <- dplyr::summarize(mean_raw_data_summary_1_wd, mean_diff = B - A, markers)
  ## Standard deviation
  sd_raw_data_summary_1 <- raw_data_lg %>%
    dplyr::group_by(i, markers)  %>%
    dplyr::summarize(standard_deviation = sd(raw_intensity))
  sd_mean_raw_data_summary_1 <- sd_raw_data_summary_1 %>%
    dplyr::group_by(markers)  %>%
    dplyr::summarize(mean_sd = mean(standard_deviation))
  ## Combine
  cb_raw_data <- dplyr::left_join(mean_diff_raw_data_summary_1_wd, sd_mean_raw_data_summary_1)
  eff_size_raw_data <- dplyr::summarize(cb_raw_data, effect_size = round(mean_diff / mean_sd, digits = 1), markers)

  return(eff_size_raw_data[,c("markers", "effect_size")])
}
