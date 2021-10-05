# Functions to generate the data

### Helper functions ###

#' Data paired.
#'
#' @details Function to check if the data is paired.
#'
#' @param data data.frame, experimental information containing (sample IDS,
#' donor IDs).
#'
#' @return logical, TRUE is the data are paired - FALSE if the data are not paired
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
#' @details Function to check that number of DE markers greater than number of
#' markers.
#'
#' @param nb_markers numeric, total number of markers.
#' @param nb_DEmarker numeric, number of differentially expressed markers.
#'
#' @return error message if the number of DE markers is greater than the total
#' number of markers.
function_DEmarkers_sup_nbmarkers <- function(nb_markers, nb_DEmarker){
  if(nb_markers < nb_DEmarker){
    stop("Number of markers is not greater than number of differentially expressed makers")
  }
}

#' Number of markers greater than 2.
#'
#' @details Function to check that the number of markers is greater than 2.
#'
#' @param nb_marker numeric, total number of markers.
#'
#' @return error message if the total number of markers is lower than 3.
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
#'
#' @return error message
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
#' @param cofactor numeric, co-factor used in the arcsinh (by default
#' cofactor = 5).
#'
#' @return data.frame of transformed data.
function_to_transform_data <- function(data, cofactor = 5) {
  asinh(data/cofactor)
}

### Function to create information ###

#' Extract marker names.
#'
#' @details Function to extract the marker names.
#'
#' @param mock_dataset data.frame, containing the cell values for each marker.
#'
#' @return vector of marker names.
function_extract_marker_names <- function(mock_dataset){
  setdiff(colnames(mock_dataset), c("group_id", "donor_id", "sample_id"))
}

#' Generate name of the DE markers.
#'
#' @details Function to generate the name of the DE markers.
#'
#' @param total_nb_marker numeric, total number of markers.
#' @param nb_DE_marker numeric, number of DE markers.
#'
#' @return vactor of marker names.
function_names_DE_markers <- function(total_nb_marker, nb_DE_marker){
  # Name
  paste("Marker", seq_len(nb_DE_marker), sep = "")
}

### Function to create data ###

#' Compute simulated cell values for one marker with markers NB informations
#'
#' @details Function to generate value for one marker with a mean and dispersion
#' specified for the negative binomiale.
#'
#' @param marker_name character, name of the marker.
#' @param nb_donor numeric, number of donors.
#' @param rho numeric, fold change.
#' @param mu0 numeric, general donor mean from which the individual mu0i will
#' be drawn.
#' @param dispersion numeric, dispersion of the markers.
#' @param subject_effect numeric, standard deviation for the normal distribution
#' from which the donor's means will be drawn (by default subject_effect = 0.01).
#' @param nb_cell_per_sample numeric, number of cells per sample (by default
#' nb_cell_per_sample = 500).
#'
#' @return data.frame of cell values.
function_value_onemarker <- function(marker_name,
                                     mu0,
                                     dispersion,
                                     subject_effect,
                                     nb_donor,
                                     nb_cell_per_sample,
                                     rho = 1){

  # Create a data.frame for experiment info
  df_experiment_info_paired <-
    data.frame(donor_id = rep(paste0("donor_", seq_len(nb_donor)), 2),
               group_id = rep(c("A", "B"), each=nb_donor),
               sample_id = paste("sample_", seq_len(nb_donor*2), sep = ""))


  # Compute mu0i for each donor
  ## Shape
  k <- mu0^2 / subject_effect^2
  ## Scale
  theta <- subject_effect^2 / mu0
  ## mu0i
  mu0i <- stats::rgamma(nb_donor, shape = k, scale = theta)
  # Compute mu_1 (mean DE marker) based on donor mean
  mu_1i <- rho*mu0i

  # Per donor
  raw_mock_dataset_expdesign <- lapply(seq_len(nb_donor), function(i_donor){
    # Cell values
    ## Condition A
    condition_A <- matrix(stats::rnbinom(n = nb_cell_per_sample, mu = mu0i[i_donor],
                                  size = dispersion))
    # Condition B
    condition_B <- matrix(stats::rnbinom(n = nb_cell_per_sample, mu = mu_1i[i_donor],
                                  size = dispersion))

    # Combine matrices
    combined_matrices <- rbind(condition_A, condition_B)
    # Raw intensity values
    combined_matrices_raw <- as.data.frame(combined_matrices)
    # Markers names
    colnames(combined_matrices_raw) <- "count" #marker_name
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
  raw_mock_dataset_expdesign <- dplyr::bind_rows(raw_mock_dataset_expdesign)
  raw_mock_dataset_expdesign$group_id <- as.factor(raw_mock_dataset_expdesign$group_id)
  raw_mock_dataset_expdesign$marker_name <- marker_name

  # Return
  return(raw_mock_dataset_expdesign)
}

#' Compute simulated cell values for one simulation with markers NB information.
#'
#' @details Function to compute the simulated cell values using a combination of
#' variable parameters, when we have prior information about the markers
#' distribution parameters (mean and dispersion of the negative binomial).
#'
#' @param variation list, list of data.frames containing the different variable
#' input paramaters to generate the data:
#'     - marker_name: name of the marker (character);
#'     - nb_donor: number of donors;
#'     - rho: fold change;
#'     - subject_effect: standard deviation for the normal distribution from
#'     which the donor's means will be drawn;
#'     - mu0: mean of the negative binomial for the gamma distribution from
#'     which the means of the different donor will be drawn;
#'     - dispersion: dispersion of the negative binomial from
#'     which the DE marker's cell values will be drawn;
#'     - nb_cell_per_sample: number of cells per sample.
#'
#' @return list with 4 slots:
#'     - df_info: data.frame of experimental information;
#'     - DEmarkers_names: vector of DE marker names;
#'     - raw_data: data.frame of raw cell values;
#'     - data: data.frame of transformed cell values.
function_create_mock_dataset_withmarkerinfo <- function(variation){
  # Generate values for all markers
  ls_raw_mock_dataset_expdesign <- lapply(variation, function(onemarker) {
    function_value_onemarker(marker_name = onemarker$marker_name,
                             mu0 = onemarker$mu0,
                             dispersion = onemarker$dispersion,
                             subject_effect = onemarker$subject_effect,
                             nb_donor = onemarker$nb_donor,
                             nb_cell_per_sample = onemarker$nb_cell_per_sample,
                             rho = onemarker$rho)
  })
  # Merge markers values
  raw_mock_dataset_expdesign_w <- dplyr::bind_rows(ls_raw_mock_dataset_expdesign)
  # Transform data
  mock_dataset <- function_to_transform_data(raw_mock_dataset_expdesign_w$count)
  mock_dataset_expdesign_w <- cbind("transformed_values" = mock_dataset,
                                  raw_mock_dataset_expdesign_w)
  # Create long formats
  raw_mock_dataset_expdesign <- raw_mock_dataset_expdesign_w %>%
    dplyr::group_by(.data$marker_name) %>%
    dplyr::mutate(rn = dplyr::row_number()) %>% # recreated unique identifier column
    tidyr::pivot_wider(names_from = .data$marker_name, values_from = .data$count)
  raw_mock_dataset_expdesign$rn <- NULL
  mock_dataset_expdesign <- mock_dataset_expdesign_w %>%
    dplyr::select(-c(.data$count)) %>%
    dplyr::group_by(.data$marker_name) %>%
    dplyr::mutate(rn = dplyr::row_number()) %>% # recreated unique identifier column
    tidyr::pivot_wider(names_from = .data$marker_name,
                       values_from = .data$transformed_values)
  mock_dataset_expdesign$rn <- NULL
  # Which markers are DE?
  df_var <- dplyr::bind_rows(variation)
  df_onlyDE <- dplyr::filter(df_var, .data$rho != 1)
  DEmarkers_names <- as.vector(df_onlyDE$marker_name)

  # Create a data.frame for experiment info
  df_experiment_info_paired <-
    unique(dplyr::select(mock_dataset_expdesign,
                         c("group_id", "donor_id", "sample_id")))

  # Create list of simulated data
  return(list("df_info" = df_experiment_info_paired,
              "DEmarkers_names" = DEmarkers_names,
              "raw_data" = raw_mock_dataset_expdesign,
              "data" = mock_dataset_expdesign))
}

### Wrapper functions ###

#' One simulation with markers NB information.
#'
#' @details Apply one simulation for a given parameter combination
#' number times (iterations), when we have prior information about the markers
#' distribution parameters (mean and dispersion of the negative binomial).
#'
#' @param variation list, list of list containing the different input
#' parameter variations @describeIn function_create_mock_dataset_withmarkerinfo.
#'
#' @return list of simulated data. Each list member contains 5 slots:
#'     - variation: the variation of input paramaters which has been used in input;
#'     - df_info: experimental information (donor IDs, group IDs, samples IDs);
#'     - DE_markers_names: name of the differentially expressed markers;
#'     - ls_mock_data: list of data.frames, each data.frame being one simulation of
#'     the cell values using the input parameters provided (list is length nb_sim).
function_apply_onesimulation_withmarkerinfo <- function(variation) {

  print(variation)
  # Apply one simulation
  output <- function_create_mock_dataset_withmarkerinfo(variation = variation)
  # Extract the experimental info which are unique for all simulation of this iteration
  df_info <- output$df_info
  # Extract mock data
  list_df_data <- output$data
  # Extract mock raw data
  list_df_raw_data <- output$raw_data
  # Extract the name(s) of the DE markers which are unique for all simulation of this iteration
  DEmarkers_names <- output$DEmarkers_names
  # Combine with the experimental info, with mock data and variation
  list_combined_output <- list("variation" = dplyr::bind_rows(variation),
                               "df_info" = df_info,
                               "DE_markers_names" = DEmarkers_names,
                               "ls_mock_data" = list_df_data,
                               "ls_mock_raw_data" = list_df_raw_data)

  # return
  return(list_combined_output)
}
