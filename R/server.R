# Server
appServer <- function(input, output, session){

  # Precomputed datasets
  observeEvent(input$p_goButton,{

    # Warning on NA
    observe({
      # Create a vector checking if all possible parameters == NA
      na_tf <- c(input$p_nb_cells == "NA",
                 input$p_nb_donors == "NA",
                 input$p_subeff == "NA",
                 input$p_mu0 == "NA",
                 input$p_dispersion == "NA")
      # Count the number of NA
      nb_na <- sum(na_tf, na.rm=TRUE)
      # Is the number of NA above 1
      nb_na_error <- nb_na > 1
      # Display message
      output$na_error <- if (nb_na_error){
        renderText({validate("ERROR: only one \"NA\" value is acceptable")})
      } else {
        NULL
      }
    })

    # Read data
    p_data <- utils::read.table(system.file("extdata",
                                            "df_precomputed_datasets.txt",
                                            package = "CyTOFpower"))

    # Plot curve
    if ( input$p_nb_donors == "NA" )  { # If no donor
      p_data_filtered <- dplyr::filter(p_data,
                                       .data$rho == input$p_fc &
                                         .data$model == input$p_modelcheckGroup &
                                         .data$nb_cell_per_sample == input$p_nb_cells &
                                         .data$subject_effect == input$p_subeff &
                                         .data$total_nb_markers == input$p_total_nb_markers &
                                         .data$mu0 == input$p_mu0 &
                                         .data$dispersion == input$p_dispersion)
      output$p_plot_pwr_nbdonors <- renderPlot({
        ggplot2::ggplot(p_data_filtered, ggplot2::aes_string(x = "nb_donor", y = "power")) +
          ggplot2::geom_point() +
          ggplot2::geom_line() +
          ggplot2::scale_y_continuous(limits = c(0, 1)) +
          ggplot2::labs(title = "Plot of power by number of donors",
                        x = "Number of donors", y = "Power") +
          ggplot2::theme_bw()
      })
    } else if ( input$p_nb_cells == "NA" ) { # If no number of cells
      p_data_filtered <- dplyr::filter(p_data,
                                       .data$rho == input$p_fc &
                                         .data$model == input$p_modelcheckGroup &
                                         .data$nb_donor == input$p_nb_donors &
                                         .data$subject_effect == input$p_subeff &
                                         .data$total_nb_markers == input$p_total_nb_markers &
                                         .data$mu0 == input$p_mu0 &
                                         .data$dispersion == input$p_dispersion)
      output$p_plot_pwr_nbcells <- renderPlot({
        ggplot2::ggplot(p_data_filtered, ggplot2::aes_string(x = "nb_cell_per_sample", y = "power")) +
          ggplot2::geom_point() +
          ggplot2::geom_line() +
          ggplot2::scale_y_continuous(limits = c(0, 1)) +
          ggplot2::labs(title = "Plot of power by number of cells",
                        x = "Number of cells per sample", y = "Power") +
          ggplot2::theme_bw()
      })
    } else if ( input$p_subeff == "NA" )  { # If no subject effect
      p_data_filtered <- dplyr::filter(p_data,
                                       .data$rho == input$p_fc &
                                         .data$model == input$p_modelcheckGroup &
                                         .data$nb_cell_per_sample == input$p_nb_cells &
                                         .data$nb_donor == input$p_nb_donors &
                                         .data$total_nb_markers == input$p_total_nb_markers &
                                         .data$mu0 == input$p_mu0 &
                                         .data$dispersion == input$p_dispersion)
      output$p_plot_pwr_subject_effect <- renderPlot({
        ggplot2::ggplot(p_data_filtered, ggplot2::aes_string(x = "subject_effect", y = "power")) +
          ggplot2::geom_point() +
          ggplot2::geom_line() +
          ggplot2::scale_y_continuous(limits = c(0, 1)) +
          ggplot2::labs(title = "Plot of power by subject effect",
                        x = "Subject Effect", y = "Power") +
          ggplot2::theme_bw()
      })
    } else if ( input$p_total_nb_markers == "NA" )  { # If no total number of markers
      p_data_filtered <- dplyr::filter(p_data,
                                       .data$rho == input$p_fc &
                                         .data$model == input$p_modelcheckGroup &
                                         .data$nb_cell_per_sample == input$p_nb_cells &
                                         .data$subject_effect == input$p_subeff &
                                         .data$nb_donor == input$p_nb_donors &
                                         .data$mu0 == input$p_mu0 &
                                         .data$dispersion == input$p_dispersion)
      output$p_plot_pwr_total_nb_markers <- renderPlot({
        ggplot2::ggplot(p_data_filtered, ggplot2::aes_string(x = "total_nb_markers", y = "power")) +
          ggplot2::geom_point() +
          ggplot2::geom_line() +
          ggplot2::scale_y_continuous(limits = c(0, 1)) +
          ggplot2::labs(title = "Plot of power by the total number of markers",
                        x = "Total number of markers", y = "Power") +
          ggplot2::theme_bw()
      })
    } else if ( input$p_mu0 == "NA" )  { # If no mu0
      p_data_filtered <- dplyr::filter(p_data,
                                       .data$rho == input$p_fc &
                                         .data$model == input$p_modelcheckGroup &
                                         .data$nb_cell_per_sample == input$p_nb_cells &
                                         .data$subject_effect == input$p_subeff &
                                         .data$nb_donor == input$p_nb_donors &
                                         .data$total_nb_markers == input$p_total_nb_markers &
                                         .data$dispersion == input$p_dispersion)
      output$p_plot_pwr_mu0 <- renderPlot({
        ggplot2::ggplot(p_data_filtered, ggplot2::aes_string(x = "mu0", y = "power")) +
          ggplot2::geom_point() +
          ggplot2::geom_line() +
          ggplot2::scale_y_continuous(limits = c(0, 1)) +
          ggplot2::labs(title = "Plot of power by NB marker mean",
                        x = "Mean mu_0", y = "Power") +
          ggplot2::theme_bw()
      })
    } else if ( input$p_dispersion == "NA" )  { # If no dispersion
      p_data_filtered <- dplyr::filter(p_data,
                                       .data$rho == input$p_fc &
                                         .data$model == input$p_modelcheckGroup &
                                         .data$nb_cell_per_sample == input$p_nb_cells &
                                         .data$subject_effect == input$p_subeff &
                                         .data$nb_donor == input$p_nb_donors &
                                         .data$total_nb_markers == input$p_total_nb_markers &
                                         .data$mu0 == input$p_mu0)
      output$p_plot_pwr_dispersion <- renderPlot({
        ggplot2::ggplot(p_data_filtered, ggplot2::aes_string(x = "dispersion", y = "power")) +
          ggplot2::geom_point() +
          ggplot2::geom_line() +
          ggplot2::scale_y_continuous(limits = c(0, 1)) +
          ggplot2::labs(title = "Plot of power by NB marker dispersion",
                        x = "Dispersion", y = "Power") +
          ggplot2::theme_bw()
      })
    } else {
      # Display value
      p_data_filtered <- dplyr::filter(p_data,
                                       .data$rho == input$p_fc &
                                         .data$model == input$p_modelcheckGroup &
                                         .data$nb_cell_per_sample == input$p_nb_cells &
                                         .data$subject_effect == input$p_subeff &
                                         .data$nb_donor == input$p_nb_donors &
                                         .data$total_nb_markers == input$p_total_nb_markers &
                                         .data$mu0 == input$p_mu0 &
                                         .data$dispersion == input$p_dispersion)

      output$p_power <- renderText(
        paste0("Power = ", dplyr::pull(p_data_filtered, .data$power)))
    }

  })

  # Personalized datasets
  # Warning on number of donors
  nb_donors_warning <- reactive({
    n <- input$nb_donors > 2
    shinyFeedback::feedbackWarning("nb_donors", !n, "At least 3 donors")
  })
  output$nb_donors_warning <- renderText(nb_donors_warning())

  # Matrix of parameters
  output$matrix <- renderUI({
    shinyMatrix::matrixInput(inputId = "param",
                             label = "Parameters per marker",
                             rows = list(
                               names = TRUE,
                               editableNames = TRUE),
                             cols = list(
                               names = TRUE,
                               editableNames = FALSE),
                             value = matrix(rep(c(1, 12, 3), input$nb_tot_markers),
                                            ncol = 3,
                                            byrow = TRUE,
                                            dimnames = list(c(paste0("m_", seq_len(input$nb_tot_markers))),
                                                            c("Fold Change", "Mean", "Dispersion")))
    )
  })

  # Show notification for no DE
  id <- NULL
  observeEvent(input$param,
               {
                 removeNotification(id)
                 fc <- any(input$param[,1] != "1")
                 if (!fc) {
                   id <<- showNotification("Please provide at least one DE marker (FC != 1)",
                                           duration = NULL,
                                           type = "error")
                 } else {
                   removeNotification(id)
                   id <<- NULL
                 }
               })

  # Once the goButton is activated
  observeEvent(input$goButton,
               {
                 # Extract matrix of parameters informations
                 df_param <- tibble::rownames_to_column(as.data.frame(input$param),
                                                        var = "marker_name")
                 colnames(df_param) <- c("marker_name", "rho", "mu0", "disp")

                 # Display a summary of the selected parameters
                 output$title_selected_var <- renderText({
                   paste("Selected parameters:", sep = "\n")
                 })
                 # Selected parameters description
                 output$selected_var <- renderUI({
                   # input$goButton
                   isolate({
                     # Introduction sentence
                     str_intro <- paste("For this simulation, you have selected
                                        the following parameters:", sep = "\n")
                     # Number of donors
                     str_nbdonor <- paste("Number of donors: ", input$nb_donors)
                     # Subject effect
                     str_subeff <- paste("Subject effect: ", input$subeff)
                     # Number of cells per sample
                     str_nbcells <- paste("Number of cells per sample: ",
                                          input$nb_cells)
                     # Total number of markers
                     str_totnbmarker <- paste("Total number of markers: ", dim(df_param)[1])
                     # Fold change
                     str_rho <- paste("Fold change of the marker(s): ", paste0(df_param$rho, collapse = "; "))
                     # Mu0
                     str_mu0 <- paste("Mean of the marker(s): ", paste0(df_param$mu0, collapse = "; "))
                     # Dispersion
                     str_disp <- paste("Dispersion of the marker(s): ", paste0(df_param$disp, collapse = "; "))
                     # Model to run
                     col_modelcheckGroup <- paste(input$modelcheckGroup, collapse = ", ")
                     str_models <- paste("Model(s) to run: ", col_modelcheckGroup)
                     # Number of simulation
                     str_sim <- paste("Number of simulation (s): ", input$nb_sim)

                     HTML(paste(str_intro, str_nbdonor,
                                str_subeff, str_nbcells,
                                str_totnbmarker,
                                # str_nbDEmarker,
                                str_rho, str_mu0, str_disp, str_models,
                                str_sim,
                                sep = '<br/>'))
                   })
                 })

                 # Generate the data
                 withProgress(message = 'Generating data', value = 0, {
                   ls_variation <- lapply(seq_len(input$nb_tot_markers), function(i){
                     data.frame(
                       "marker_name" = df_param[i, "marker_name"],
                       "mu0" = as.numeric(df_param[i, "mu0"]),
                       "dispersion" = as.numeric(df_param[i, "disp"]),
                       "subject_effect" = input$subeff,
                       "nb_donor" = input$nb_donors,
                       "nb_cell_per_sample" = input$nb_cells,
                       # "rho" = as.numeric(df_rho()[i]))
                       "rho" = as.numeric(df_param[i, "rho"]))
                   })
                   # Pipeline
                   ls_res <- lapply(seq_len(input$nb_sim), function(sim_id){
                     # (1) generate a single dataset
                     sim_data <- function_apply_onesimulation_withmarkerinfo(ls_variation)
                     # (2) apply methods
                     res_mod <- function_apply_modelcomputations_modelchoice(list_combined_output = sim_data,
                                                                             model = input$modelcheckGroup)
                     # (3) store results and discard data
                     # Observed variance/effect size
                     # Merge counts data
                     raw_data <- cbind(sim_data$ls_mock_raw_data, "i" = sim_id)
                     # Long format
                     raw_data_lg <- tidyr::pivot_longer(raw_data,
                                                        # cols = starts_with("Marker"),
                                                        cols = -c("i", "group_id", "donor_id", "sample_id"),
                                                        names_to = "markers",
                                                        values_to = "raw_intensity")
                     # Variance and effect size
                     obs_variance <- compute_variance(raw_data_lg)
                     obs_effectsize <- compute_effectsize(raw_data_lg)
                     # Discard
                     rm(sim_data)
                     # Return
                     return(list("res_models" = res_mod,
                                 "obs_variance" = obs_variance,
                                 "obs_effectsize" = obs_effectsize))
                     # (4) repeat
                   })
                 })

                 # Combine results
                 # Models
                 df_res_models <- do.call("bind_rows", list(lapply(ls_res, function(x){
                   x[["res_models"]]
                 }), .id = "i"))
                 # Variance
                 df_obs_variance_persim <- do.call("bind_rows", list(lapply(ls_res, function(x){
                   x[["obs_variance"]]
                 }), .id = "i"))
                 ## mean
                 df_obs_variance <- df_obs_variance_persim %>%
                   dplyr::group_by(.data$donor_id) %>%
                   dplyr::summarize(variance = round(mean(.data$variance), digits = 2))
                 # Effect size
                 df_obs_effectsize_persim <- do.call("bind_rows", list(lapply(ls_res, function(x){
                   x[["obs_effectsize"]]
                 }), .id = "i"))
                 ## Mean
                 df_obs_effectsize <- df_obs_effectsize_persim %>%
                   dplyr::group_by(.data$markers) %>%
                   dplyr::summarize(effect_size = round(mean(.data$effect_size), digits = 1),
                                    observed_FC = round(mean(.data$observed_FC), digits = 1))

                 # Name of DE(s) markers
                 df_var <- dplyr::bind_rows(ls_variation)
                 df_var_DE <- dplyr::filter(df_var, .data$rho != 1)
                 output$DE_names <- renderText(
                   paste("Name(s) of the DE markers:", paste(df_var_DE$marker_name, collapse = "; "))
                 )


                 # Output
                 output$title_var <- renderText(paste("Observed sample variance (mean across simulations)"))
                 output$tab_obs_variance <- DT::renderDataTable({df_obs_variance}, options = list(pageLength = 5))
                 output$title_effsize <- renderText(paste("Observed Cohen's effect size"))
                 output$tab_obs_effectsize <- DT::renderDataTable({df_obs_effectsize}, options = list(pageLength = 5))
                 output$results_models <- DT::renderDataTable({df_res_models})

                 # Compute power
                 withProgress(message = "Compute power", value = 0, {
                   pwr_values <- compute_pwr(df_res_models)
                 })
                 colnames(pwr_values) <- c("model", "markers", "power")
                 output$title_power <- renderText(paste("Power"))
                 output$results_pwr_values <- DT::renderDataTable({pwr_values[,c("model", "markers", "power")]})

                 # Effect size - power
                 ## Combine
                 pwr_eff <- dplyr::left_join(pwr_values, df_obs_effectsize, by = "markers")
                 ## Plot
                 if(length(input$modelcheckGroup) > 1){
                   output$plot_pwr_eff <- renderPlot({
                     ggplot2::ggplot(pwr_eff, ggplot2::aes_string(x = "effect_size", y = "power", col = "model")) +
                       ggplot2::geom_point(size = 2,
                                           position = ggplot2::position_jitter(width = 0.01, height = 0.01)) + #position = "jitter" position = position_jitter()
                       ggplot2::ylim(0, 1.3) +
                       ggplot2::labs(title="Plot of power by effect size",
                                     x ="Cohen effect size", y = "Power") +
                       ggplot2::theme_bw()
                   })
                 } else{
                   output$plot_pwr_eff <- renderPlot({
                     ggplot2::ggplot(pwr_eff, ggplot2::aes_string(x = "effect_size", y = "power", col = "model")) +
                       ggplot2::geom_point(size = 2) + #position = "jitter" position = position_jitter()
                       ggplot2::ylim(0, 1) +
                       ggplot2::labs(title="Plot of power by effect size",
                                     x ="Cohen effect size", y = "Power") +
                       ggplot2::theme_bw()
                   })
                 }

                 # Wide format
                 pwr_eff_w <- tidyr::pivot_wider(data = pwr_eff,
                                                 id_cols = c("markers", "effect_size", "observed_FC"),
                                                 names_from = "model",
                                                 values_from = "power")
                 ## Add non-DE marker
                 nonDE_obs_effectsize <- df_obs_effectsize[!(df_obs_effectsize$markers %in% pwr_eff_w$markers),]
                 col_to_add <- colnames(pwr_eff_w)[!(colnames(pwr_eff_w) %in% colnames(nonDE_obs_effectsize))]
                 nonDE_obs_effectsize[,col_to_add] <- NA
                 pwr_eff_w_nonDE <- rbind(pwr_eff_w, nonDE_obs_effectsize)
                 output$results_pwr_with_effsize <- DT::renderDataTable(
                   # DT::datatable({pwr_eff_w_nonDE},
                   pwr_eff_w_nonDE,
                   # class = 'cell-border stripe',
                   # extensions = 'Buttons',
                   # options = list(autoWidth = TRUE,
                   #                dom = 'Bfrtip',
                   #                buttons = c("cvs"))),
                   server = FALSE)
               }
  )

}
