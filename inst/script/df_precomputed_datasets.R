# Generation of the df_precomputed_datasets.txt file
#
# This file was generated using un-exported functions of the CyTOFpower package.
# Each set of parameter was generated as described below:

# A function regrouping the different steps of the pipeline was created:
fct_workflow <- function(nb_sim, ls_variation, modelcheckGroup){
  lapply(1:nb_sim, function(sim_id){
    # (1) Generate a single dataset
    sim_data <- CyTOF:::function_apply_onesimulation_withmarkerinfo(ls_variation)
    # (2) Apply methods
    res_mod <- CyTOF:::function_apply_modelcomputations_modelchoice(
      list_combined_output = sim_data,
      model = modelcheckGroup)
    # (3) Discard data for memory space
    rm(sim_data)
    # Return
    return(list("res_models" = res_mod))
    # (4) Repeat
  })
}

# Each parameter for one set of parameters was defined:
# Number of simulations
nb_sim
# Number of donors
nb_donors
# Number of markers
nb_markers
## Parameters of markers
name_marker <- paste0("m_", 1:nb_markers)
## Log Change
vect_logchange
logchange <- c(vect_logchange, rep(1, nb_markers-length(vect_logchange)))
## Subject effect
value_subeff
subeff <- rep(value_subeff, nb_markers)
## Mean
mean_negbino
## Dispersion
disp_negbino
# Models to test
modelcheck
# Number of cells per sample
nb_cells
# Create a data.frame containing these parameters
df_var <- data.frame(marker_name = name_marker,
                     mu0 = rep(mean_negbino, nb_markers),
                     dispersion = rep(disp_negbino, nb_markers),
                     subject_effect = subeff,
                     nb_donor = rep(nb_donors, nb_markers),
                     nb_cell_per_sample = rep(nb_cells, nb_markers),
                     rho = logchange)
# Convert the parameters data.frame as a list
ls_var <- split(df_var, df_var$marker_name)

# The workflow was applied
# Workflow
ls_res <- fct_workflow(nb_sim = nb_sim,
                       ls_variation = ls_var,
                       modelcheckGroup = modelcheck)

# Results for the different simulations were combined
# Models
df_res_models <- do.call("bind_rows", list(lapply(ls_res, function(x){
  x[["res_models"]]
}), .id = "i"))
# Power was computed for alpha = 0.05
pwr_values <- CyTOFpower::compute_pwr(df_res_models)
colnames(pwr_values)[2] <- "marker_name"
pwr_values <- left_join(pwr_values,
                        df_var)

