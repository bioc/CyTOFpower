# Shiny

# UI
myui <- fluidPage(
  # App title
  titlePanel("Power simulation"),

  # Sidebar panel
  sidebarLayout(

    # Sidebar panel for inputs ----
    sidebarPanel(

      # Input: Field to enter the number of donors ----
      numericInput(inputId = "nb_donors",
                   label = "Number of donors:",
                   value = 10),

      # Input: Slider for the subject effect ----
      sliderInput(inputId = "subeff",
                  label = "Subject effect:",
                  min = 0.01,
                  max = 3,
                  value = 0.5),

      # Input: Field to enter the number of cells per sample ----
      numericInput(inputId = "nb_cells",
                   label = "Number of cells per samples:",
                   value = 500),

      # Input: Slider for the number of markers ----
      sliderInput(inputId = "nb_tot_markers",
                  label = "Total number of markers:",
                  min = 3,
                  max = 50,
                  value = 30),

      # Input: Slider window for the number of DE markers ----
      sliderInput(inputId = "nb_DEmarkers",
                  label = "Number of DE markers:",
                  min = 1,
                  max = 50,
                  value = 1),

      # Input: Slider for the fold change ----
      # sliderInput(inputId = "rho",
      #             label = "Fold change of the DE markers:",
      #             min = 0.01,
      #             max = 3,
      #             value = 1.2),
      uiOutput("sliders_nb_DEmarkers"),

      # Input: Check box to choose which model to run ----
      checkboxGroupInput(inputId = "modelcheckGroup",
                         label = "Models to run",
                         choices = list("CytoGLMM - GLMM" = "cytoglmm",
                                        "CytoGLMM - GLM with bootstrap" = "cytoglm",
                                        "diffcyt - limma with fixed effect" = "testDS_limma_fixed",
                                        "diffcyt - limma with random effect" = "testDS_limma_random",
                                        "diffcyt - lme4" = "testDS_lmm"),
                         selected = "cytoglmm"),

      # Input: Field to enter the number of simulations ----
      numericInput(inputId = "nb_sim",
                   label = "Number of simulations for the set of selected parameters:",
                   value = 10,
                   min = 1,
                   max = 100),

      # Input: Action button to submit the paramaters
      actionButton(inputId = "goButton",
                   label = "Run simulation")
    ),
    # Main panel
    mainPanel("The main panel will return the power of the simulation, with the observed
                variance as well as the effect size",
              # h4("Selected fold change(s):"),
              # tableOutput(outputId = "table_rho_m"),
              #textOutput(outputId = "toprint_rho_m"),
              # uiOutput: selected parameters ----
              h4(textOutput(outputId = "title_selected_var")),
              uiOutput(outputId = "selected_var"),
              # Dataset
              h4(textOutput(outputId = "simulated_data")),
              #DT::dataTableOutput(outputId = "df_info"),
              textOutput(outputId = "DE_names"),
              # Power
              h4(textOutput(outputId = "title_power")),
              DT::dataTableOutput(outputId = "results_pwr_values"),
              # Observed variance and effect size
              h4(textOutput(outputId = "title_var")),
              DT::dataTableOutput(outputId = "tab_obs_variance"),
              h4(textOutput(outputId = "title_effsize")),
              DT::dataTableOutput(outputId = "tab_obs_effectsize"),
              ## Results
              DT::dataTableOutput(outputId = "results_models")


    )

  )
)


# Server
myserver = function(input, output, session){

  # Adapt maximum of the number of DE marker which cannot be higher than the
  # total number of markers
  observe({
    nb_markers <- input$nb_tot_markers
    updateSliderInput(session,
                      inputId = "nb_DEmarkers",
                      max = input$nb_tot_markers)
  })

  # Create the fold change per marker
  output$sliders_nb_DEmarkers <- renderUI({
    nb_DEmarkers <- input$nb_DEmarkers
    lapply(seq_len(nb_DEmarkers), function(i) {
      sliderInput(inputId = paste0("rho_m", i),
                  label = paste0("Fold change of DE marker ", i),
                  min = 0.01,
                  max = 3,
                  value = 1.2)
    })
  })

  # Print fold change for the different marker(s)
  ## TODO: integrate in the observeEvent gobutton
  ## or put in the sidebarLayout https://stackoverflow.com/questions/46958390/is-it-possible-to-put-a-data-table-in-the-sidebar-in-shinydashboard
  observe({
    df_rho <<- reactive({
      # Get values of the FC
      ls <- lapply(seq_len(input$nb_DEmarkers), function(i) {
        input[[paste0("rho_m", i)]]
      })
      names(ls) <- paste0("Marker", seq_len(input$nb_DEmarkers))
      bind_rows(ls, .id = "marker")
    })
    # output$table_rho_m <- renderTable({df_rho()})
    # output$toprint_rho_m <- renderText({paste0(df_rho(), collapse = "; ")})
  })

  # Once the goButton is activated
  observeEvent(input$goButton,
               {
                 output$title_selected_var <- renderText({
                   paste("Selected parameters:", sep = "\n")
                 })
                 # Selected parameters description
                 output$selected_var <- renderUI({
                   # input$goButton
                   isolate({
                     # Introduction sentence
                     str_intro <- paste("For this simulation, you have selected the following parameters:", sep = "\n")
                     # Number of donors
                     str_nbdonor <- paste("Number of donors: ", input$nb_donors)
                     # Subject effect
                     str_subeff <- paste("Subject effect: ", input$subeff)
                     # Number of cells per sample
                     str_nbcells <- paste("Number of cells per sample: ", input$nb_cells)
                     # Total number of markers
                     str_totnbmarker <- paste("Total number of markers: ", input$nb_tot_markers)
                     # Number of DE markers
                     str_nbDEmarker <- paste("Number of DE markers: ", input$nb_DEmarkers)
                     # Fold change
                     str_rho <- paste("Fold change for the DE marker(s): ", paste0(df_rho(), collapse = "; ")) #input$rho)
                     # Model to run
                     col_modelcheckGroup <- paste(input$modelcheckGroup, collapse = ", ")
                     str_models <- paste("Model(s) to run: ", col_modelcheckGroup)
                     # Number of simulation
                     str_sim <- paste("Number of simulation (s): ", input$nb_sim)

                     HTML(paste(str_intro, str_nbdonor,
                                str_subeff, str_nbcells,
                                str_totnbmarker, str_nbDEmarker,
                                str_rho, str_models,
                                str_sim,
                                sep = '<br/>'))
                   })
                 })

                 # Generate the data
                 withProgress(message = 'Generating data', value = 0, {
                   ls_variation <- list(nb_DE_marker = input$nb_DEmarkers,
                                        rho = as.vector(t(df_rho())), #c(1.2), #input$rho,
                                        nb_donor = input$nb_donors,
                                        subject_effect = input$subeff,
                                        total_nb_marker = input$nb_tot_markers,
                                        nb_cell_per_sample = input$nb_cells)
                   sim_data <- function_wrapper_apply_simulation_nbtimes_multiFC(variation = ls_variation,
                                                                                 nb_sim =  input$nb_sim)
                 })
                 output$simulated_data <- renderText(
                   paste("Dataset generated:", length(sim_data))
                 )
                 #output$df_info <- DT::renderDataTable({sim_data[[1]]$df_info})
                 output$DE_names <- renderText(
                   paste("Name(s) of the DE markers:", paste(sim_data$DE_markers_names, collapse = "; "))
                 )

                 # Merge counts data
                 raw_data <- bind_rows(sim_data$ls_mock_raw_data, .id = "i")
                 # Long format
                 raw_data_lg <- tidyr::pivot_longer(raw_data,
                                                    cols = starts_with("Marker"),
                                                    names_to = "markers",
                                                    values_to = "raw_intensity")
                 # Variance and effect size
                 obs_variance <- compute_variance(raw_data_lg)
                 obs_effectsize <- compute_effectsize(raw_data_lg)
                 output$title_var <- renderText(paste("Observed sample variance"))
                 output$tab_obs_variance <- DT::renderDataTable({obs_variance}, options = list(pageLength = 5))
                 output$title_effsize <- renderText(paste("Cohen's effect size"))
                 output$tab_obs_effectsize <- DT::renderDataTable({obs_effectsize}, options = list(pageLength = 5))

                 # Run the models
                 withProgress(message = "Running the model(s)", value = 0, {
                   res_mod <- function_wrapper_apply_modelcomputations_nbtimes_modelchoice(list_combined_output = sim_data,
                                                                                           model = input$modelcheckGroup)
                 })
                 output$results_models <- DT::renderDataTable({res_mod})

                 # Compute power
                 withProgress(message = "Compute power", value = 0, {
                   pwr_values <- compute_pwr(res_mod)
                 })
                 output$title_power <- renderText(paste("Power"))
                 output$results_pwr_values <- DT::renderDataTable({pwr_values[,c("model", "marker_id", "power")]})
               }
  )


}

# App
app <- shinyApp(ui = myui, server = myserver)