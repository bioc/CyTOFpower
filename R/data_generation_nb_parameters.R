#' Compute simulated cell values for one marker with markers NB informations
#'
#' @details Function to generate value for one marker with a mean and dispersion specified
#' for the negative binomiale.
#'
#' @param marker_name character, name of the marker.
#' @param nb_donor numeric, number of donors.
#' @param rho numeric, fold change.
#' @param mu0 numeric, general donor mean from which the individual mu0i will
#' be drawn.
#' @param dispersion_marker numeric, dispersion of the markers.
#' @param subject_effect numeric, standard deviation for the normal distribution from
#' which the donor's means will be drawn (by default subject_effect = 0.01).
#' @param nb_cell_per_sample numeric, number of cells per sample (by default nb_cell_per_sample = 500).
function_value_onemarker <- function(marker_name,
                                     mu0,
                                     dispersion,
                                     subject_effect,
                                     nb_donor,
                                     nb_cell_per_sample,
                                     rho = 1){

  # Create a data.frame for experiment info
  df_experiment_info_paired <-
    data.frame(donor_id = rep(paste0("donor_", 1:nb_donor), 2),
               group_id = rep(c("A", "B"), each=nb_donor),
               sample_id = paste("sample_", 1:(nb_donor*2), sep = ""))


  # Compute mu0i for each donor
  ## Shape
  k <- mu0^2 / subject_effect^2
  ## Scale
  theta <- subject_effect^2 / mu0
  ## mu0i
  mu0i <- rgamma(nb_donor, shape = k, scale = theta)
  # Compute mu_1 (mean DE marker) based on donor mean
  mu_1i <- rho*mu0i

  # Per donor
  raw_mock_dataset_expdesign <- plyr::adply(1:nb_donor, 1, function(i_donor){
    # Cell values
    ## Condition A
    condition_A <- matrix(rnbinom(n = nb_cell_per_sample, mu = mu0i[i_donor],
                                  size = dispersion))
    # Condition B
    condition_B <- matrix(rnbinom(n = nb_cell_per_sample, mu = mu_1i[i_donor],
                                  size = dispersion))

    # Combine matrices
    combined_matrices <- rbind(condition_A, condition_B)
    # Raw intensity values
    combined_matrices_raw <- as.data.frame(combined_matrices)
    # Markers names
    colnames(combined_matrices_raw) <- marker_name
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

  # Return
  return(raw_mock_dataset_expdesign)
}

# Small dataset creation
ls_var <- list(data.frame(marker_name = "test1",
                 mu0 = 10,
                 dispersion = 0.5,
                 subject_effect = 0.3,
                 nb_donor = 2,
                 nb_cell_per_sample = 20,
                 rho = 1),
               data.frame(marker_name = "test2",
                 mu0 = 10,
                 dispersion = 0.5,
                 subject_effect = 0.3,
                 nb_donor = 2,
                 nb_cell_per_sample = 20,
                 rho = 3))

#' Compute simulated cell values for one simulation with markers NB information.
#'
#' @details Function to compute the simulated cell values using a combination of
#' variable parameters, when we have prior information about the markers distribution
#' parameters (mean and dispersion of the negative binomial).
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
  raw_mock_dataset_expdesign <- plyr::join_all(ls_raw_mock_dataset_expdesign,
                                               by = c("group_id", "donor_id", "sample_id"))
  # Transform data
  mock_dataset <- function_to_transform_data(dplyr::select(raw_mock_dataset_expdesign,
                                                           -c("group_id", "donor_id", "sample_id")))
  mock_dataset_expdesign <- cbind(mock_dataset,
                                  dplyr::select(raw_mock_dataset_expdesign,
                                                c("group_id", "donor_id", "sample_id")))
  # Which markers are DE?
  ls_isDE <- lapply(variation, function(onemarker){
    # DE marker if rho != 1
    if (onemarker$rho != 1) {
      DEmarker <-  TRUE
    } else {
      DEmarker <- FALSE
    }
    # Return
    return(data.frame("marker_name" = onemarker$marker_name,
                      "isDEmarker" = DEmarker))
  })
  df_isDE <- bind_rows(ls_isDE)
  ## Select DE markers
  df_onlyDE <- dplyr::filter(df_isDE, isDEmarker == TRUE)
  DEmarkers_names <- as.vector(df_onlyDE$marker_name)
  # Create list of simulated data
  return(list("df_info" = df_experiment_info_paired,
              "DEmarkers_names" = DEmarkers_names,
              "raw_data" = raw_mock_dataset_expdesign,
              "data" = mock_dataset_expdesign))
}

### Wrapper functions ###

#' Wrapper to repeat simulations with markers NB information.
#'
#' @details Function to repeat simulations for a given parameter combination
#' number times (iterations), when we have prior information about the markers distribution
#' parameters (mean and dispersion of the negative binomial).
#'
#' @param variation list, list of list containing the different input
#' parameter variations @describeIn function_create_mock_dataset_withmarkerinfo.
#' @param nb_sim numeric, number of simulation.
#' @return list of simulated data. Each list member contains 5 slots:
#'     - variation: the variation of input paramaters which has been used in input;
#'     - df_info: experimental information (donor IDs, group IDs, samples IDs);
#'     - DE_markers_names: name of the differentially expressed markers;
#'     - ls_mock_data: list of data.frames, each data.frame being one simulation of
#'     the cell values using the input parameters provided (list is length nb_sim).
function_wrapper_apply_simulation_nbtimes_withmarkerinfo <- function(variation, nb_sim = 10) {

  print(variation)
  # Repeat simulation nb_sim times
  output <- replicate(nb_sim, function_create_mock_dataset_withmarkerinfo(variation = variation),
                      simplify = FALSE)
  # Extract the experimental info which are unique for all simulation of this iteration
  df_info <- (unique(lapply(output, function(x) x$df_info)))[[1]]
  # Extract mock data and add the iteration nb
  list_df_data <- lapply(output, function(x) x$data)
  # Extract mock raw data and add the iteration nb
  list_df_raw_data <- lapply(output, function(x) x$raw_data)
  # Extract the name(s) of the DE markers which are unique for all simulation of this iteration
  DEmarkers_names <- (unique(lapply(output, function(x) x$DEmarkers_names)))[[1]]
  # Combine with the experimental info, with mock data and variation
  list_combined_output <- list("variation" = bind_rows(variation),
                               "df_info" = df_info,
                               "DE_markers_names" = DEmarkers_names,
                               "ls_mock_data" = list_df_data,
                               "ls_mock_raw_data" = list_df_raw_data)

  # return
  return(list_combined_output)
}
