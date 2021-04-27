# Functions to generate the data

### Helper functions ###

#' Data paired.
#'
#' @details Function to check if the data is paired.
#'
#' @param data data.frame, experimental information containing (sample IDS,
#' donor IDs).
function_is_data_paired <- function(data){
  # Number of samples
  nb_samples <- length(unique(data$sample_id))
  # Number of donors
  nb_donors <- length(unique(data$donor_id))
  # Are the data paired?
  if(nb_samples == nb_donors){
    message("The dataset is not paired: please use bootstrap GLMM model")
    are_paired_check <- "FALSE"
  } else {
    message("Your dataset is paired, you can used the GLMM model or diffcyt package")
    are_paired_check <- "TRUE"
  }

  return(are_paired_check)
}


#' Check on the number of DE markers.
#'
#' @details Function to check that number of DE markers greater than number of markers.
#'
#' @param nb_markers numeric, total number of markers.
#' @param nb_DEmarker numeric, number of differentially expressed markers.
function_DEmarkers_sup_nbmarkers <- function(nb_markers, nb_DEmarker){
  if(nb_markers < nb_DEmarker){
    stop("Number of markers is not greater than number of differentially expressed makers")
  }
}

#' Number of markers greater than 2.
#'
#' @details Function to check that the number of markers is greater than 2
#'
#' @param nb_marker numeric, total number of markers.
function_check_nbmarkers <- function(nb_marker){
  if(nb_marker < 3){
    stop("Number of markers should be greater than 2")
  }
}

#' Number of DE markers greater than 1.
#'
#' @details Function to check that the number of DE markers is greater than 1.
#'
#' @param nb_DEmarker numeric, number of differentially expressed markers.
function_check_nbDEmarkers <- function(nb_DEmarker){
  if(nb_DEmarker < 1){
    stop("Number of markers should be greater than 0")
  }
}

#' arcsinh transformation.
#'
#' @details Function to transform the data with the recommended transformation
#' for cyTOF data: arcsinh with cofactor equals to 5.
#'
#' @param data data.frame, cell values to transform.
#' @param cofactor numeric, co-factor used in the arcsinh (by default cofactor = 5).
function_to_transform_data <- function(data, cofactor = 5) {
  asinh(data/cofactor)
}

### Function to create information ###

#' Extract marker names.
#'
#' @details Function to extract the marker names.
#'
#' @param mock_dataset data.frame, containing the cell values for each marker.
function_extract_marker_names <- function(mock_dataset){
  # Get the index
  idx <- grep("Marker", colnames(mock_dataset))
  return(colnames(mock_dataset)[idx])
}

#' Generate name of the DE markers.
#'
#' @details Function to generate the name of the DE markers.
#'
#' @param total_nb_marker numeric, total number of markers.
#' @param nb_DE_marker numeric, number of DE markers.
function_names_DE_markers <- function(total_nb_marker, nb_DE_marker){
  # Number
  #nb <- seq(from = (total_nb_marker - nb_DE_marker)+1, to = total_nb_marker)
  # Name
  paste("Marker", 1:nb_DE_marker, sep = "")
}

### Function to create data ###

#' Create mock dataset.
#'
#' @details Function to create a mock dataset with subject effect in addition with
#' the possibility to vary the number of donors and rho.
#'
#' @param nb_donor numeric, number of donors.
#' @param rho numeric, fold change.
#' @param mu_0 numeric, general donor mean from which the individual mu_0i will
#' be drawn (by default mu_0 = 0).
#' @param dispersion_marker numeric, dispersion of the markers
#' (by default dispersion_marker = 0.5).
#' @param subject_effect numeric, standard deviation for the normal distribution from
#' which the donor's means will be drawn (by default subject_effect = 0.01).
#' @param marker_effect_size numeric, standard deviation for the normal distribution
#' from which the DE marker's mean will be drawn (by default marker_effect_size = 0).
#' @param nb_markers numeric, total number of markers (by default nb_markers = 30).
#' @param nb_DE_marker numeric, number of differentially expressed markers
#' (by default nb_DE_marker = 1).
#' @param nb_cell_per_sample numeric, number of cells per sample (by default nb_cell_per_sample = 500).
function_create_mock_dataset_withsubjeff_vary_nandmultiFC <- function(nb_donor,
                                                                      rho,
                                                                      mu_0 = 12,
                                                                      dispersion_marker = 0.5,
                                                                      subject_effect = 0.01,
                                                                      #marker_effect_size = 0,
                                                                      nb_cell_per_sample = 500,
                                                                      nb_markers = 30,
                                                                      nb_DE_marker = 1){
  # Number of marker
  ## Check number of markers greater than 2
  function_check_nbmarkers(nb_marker = nb_markers)
  ## Check number of DE markers greater than 1
  function_check_nbDEmarkers(nb_DEmarker = nb_DE_marker)
  ## Check number of DE markers greater than number of markers
  function_DEmarkers_sup_nbmarkers(nb_markers = nb_markers,
                                   nb_DEmarker = nb_DE_marker)

  # Create a data.frame for experiment info
  df_experiment_info_paired <-
    data.frame(donor_id = rep(paste0("donor_", 1:nb_donor), 2),
               group_id = rep(c("A", "B"), each=nb_donor),
               sample_id = paste("sample_", 1:(nb_donor*2), sep = ""))
  ## Shape
  k <- mu_0^2 / subject_effect^2
  ## Scale
  theta <- subject_effect^2 / mu_0
  ## mu_0i
  mu_0i <- rgamma(nb_donor, shape = k, scale = theta)
  # Compute mu_1 (mean DE marker) based on donor mean
  if( length(rho) > 1){
    mu_0i_rep <- matrix(mu_0i, ncol = nb_donor, nrow = length(rho), byrow = TRUE)
    mu_1i <- rho*mu_0i_rep
  } else {
    mu_1i <- rho*mu_0i
  }

  # Per donor
  raw_mock_dataset_expdesign <- plyr::adply(1:nb_donor, 1, function(i_donor){
    # Cell values
    ## Condition A
    condition_A <- matrix(rnbinom(n = nb_cell_per_sample*nb_markers, mu = mu_0i[i_donor], size = dispersion_marker),
                          ncol=nb_markers)
    # Condition B
    if( length(rho) > 1) {
      de_values <- sapply(1:length(rho), function(i_rho){
        matrix(rnbinom(n = nb_cell_per_sample, mu = mu_1i[i_rho,i_donor], size = dispersion_marker), ncol = 1)
      })
    } else {
      de_values <- matrix(rnbinom(n = nb_cell_per_sample*nb_DE_marker, mu = mu_1i[i_donor], size = dispersion_marker), ncol = nb_DE_marker)
    }
    condition_B <- cbind(
      de_values,
      matrix(rnbinom(n = nb_cell_per_sample*(nb_markers-nb_DE_marker), mu = mu_0i[i_donor], size = dispersion_marker),
             ncol=(nb_markers - nb_DE_marker))
    )

    # Combine matrices
    combined_matrices <- rbind(condition_A, condition_B)
    # Raw intensity values
    combined_matrices_raw <- as.data.frame(combined_matrices)
    # Markers names
    colnames(combined_matrices_raw) <- paste0("Marker", seq_len(nb_markers))
    # Add the conditions
    combined_matrices_raw$"group_id" <- rep(c("A", "B"), each = nb_cell_per_sample)
    # Add the donor ID
    combined_matrices_raw$"donor_id" <- unique(df_experiment_info_paired$donor_id)[i_donor]
    # Add sample ID
    combined_matrices_raw_expdesign <- dplyr::left_join(combined_matrices_raw,
                                                        df_experiment_info_paired,
                                                        by = c("donor_id", "group_id"))
    return(combined_matrices_raw_expdesign)
  })
  raw_mock_dataset_expdesign$group_id %<>% as.factor
  raw_mock_dataset_expdesign$X1 <- NULL

  # Generate the name of the DE markers
  ## I always take the last ones:
  ## i.e. from (total_nb_marker - nb_DE_marker)+1 to total_nb_marker
  DEmarkers_names <- function_names_DE_markers(total_nb_marker = nb_markers,
                                               nb_DE_marker = nb_DE_marker)
  # Transform the data
  mock_dataset <- function_to_transform_data(dplyr::select(raw_mock_dataset_expdesign,
                                                           starts_with("Marker")))
  mock_dataset_expdesign <- cbind(mock_dataset,
                                  dplyr::select(raw_mock_dataset_expdesign, !starts_with(("Marker"))))

  return(list("df_info" = df_experiment_info_paired,
              "DEmarkers_names" = DEmarkers_names,
              "raw_data" = raw_mock_dataset_expdesign,
              "data" = mock_dataset_expdesign))
}

### Wrapper functions ###

#' Compute simulation.
#'
#' @details Function to compute the simulation using a combination of variable
#' parameters.
#'
#' @param variation data.frame, data.frame containing the different variable
#' input paramaters to generate the data:
#'     - nb_donor: number of donors;
#'     - rho: fold change;
#'     - subject_effect: standard deviation for the normal distribution from
#'     which the donor's means will be drawn;
#'     - marker_effect_size: standard deviation for the normal distribution from
#'     which the DE marker's means will be drawn;
#'     - nb_DE_marker: number of differentially expressed markers;
#'     - total_nb_marker: total number of markers;
#'     - nb_cell_per_sample: number of cells per sample.
function_to_compute_simulations_multiFC <- function(variation){
  # Create mock datasets
  ls_sim <- function_create_mock_dataset_withsubjeff_vary_nandmultiFC(nb_donor = variation$nb_donor,
                                                                      rho = variation$rho,
                                                                      subject_effect = variation$subject_effect,
                                                                      #marker_effect_size = variation$marker_effect_size,
                                                                      nb_DE_marker = variation$nb_DE_marker,
                                                                      nb_markers = variation$total_nb_marker,
                                                                      nb_cell_per_sample = variation$nb_cell_per_sample)

  # Return
  return(ls_sim)
}

#' Wrapper to repeat simulations.
#'
#' @details Function to repeat simulations for a given parameter combination
#' number times (iterations).
#'
#' @param variation list, list of data.frame containing the different input.
#' parameter variations @describeIn function_to_compute_simulations_multiFC.
#' @param nb_sim numeric, number of simulation.
#' @return list of simulated data. Each list member contains 5 slots:
#'     - variation: the variation of input paramaters which has been used in input;
#'     - df_info: experimental information (donor IDs, group IDs, samples IDs);
#'     - DE_markers_names: name of the differentially expressed markers;
#'     - ls_mock_data: list of data.frames, each data.frame being one simulation of
#'     the cell values using the input parameters provided (list is length nb_sim).
function_wrapper_apply_simulation_nbtimes_multiFC <- function(variation, nb_sim = 10) {

  print(variation)
  # Repeat simulation nb_sim times
  output <- replicate(nb_sim, function_to_compute_simulations_multiFC(variation = variation), simplify = FALSE)
  # Extract the experimental info which are unique for all simulation of this iteration
  df_info <- (unique(lapply(output, function(x) x$df_info)))[[1]]
  # Extract mock data and add the iteration nb
  list_df_data <- lapply(output, function(x) x$data)
  # Extract mock raw data and add the iteration nb
  list_df_raw_data <- lapply(output, function(x) x$raw_data)
  # Extract the name(s) of the DE markers which are unique for all simulation of this iteration
  DEmarkers_names <- (unique(lapply(output, function(x) x$DEmarkers_names)))[[1]]
  # Combine with the experimental info, with mock data and variation
  list_combined_output <- list("variation" = variation,
                               "df_info" = df_info,
                               "DE_markers_names" = DEmarkers_names,
                               "ls_mock_data" = list_df_data,
                               "ls_mock_raw_data" = list_df_raw_data)

  # return
  return(list_combined_output)
}
