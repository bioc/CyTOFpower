# UI
appUI <- fluidPage(
  # For warnings
  shinyFeedback::useShinyFeedback(),
  shinyjs::useShinyjs(),

  navbarPage("CyTOFpower",

             # Precomputed datasets
             tabPanel("Precomputed dataset",
                      "This panel allows to search a grid of parameters that have been
                      precomputed. The value NA displays a power curve for the different
                      values of this parameter.",
                      # Sidebar panel
                      sidebarLayout(
                        # Sidebar panel for inputs ----
                        sidebarPanel(

                          # Input: Select the number of donors ----
                          selectInput(inputId = "p_nb_donors",
                                      label = "Number of donors:",
                                      choices = list_param_precomupteddataset$nb_donor),
                          # textOutput("na_error"),

                          # Input: Select the subject effect ----
                          selectInput(inputId = "p_subeff",
                                      label = "Subject effect:",
                                      choices = list_param_precomupteddataset$subject_effect),

                          # Input: Select the number of cells per sample ----
                          selectInput(inputId = "p_nb_cells",
                                      label = "Number of cells per sample:",
                                      choices = list_param_precomupteddataset$nb_cell_per_sample),

                          # Input: Select the total number of markers ----
                          selectInput(inputId = "p_total_nb_markers",
                                      label = "Total number of markers:",
                                      choices = list_param_precomupteddataset$total_nb_markers),

                          # Input: Select the fold change
                          selectInput(inputId = "p_fc",
                                      label = "Marker's Fold Change (FC)",
                                      choice = c(1.1, 1.2, 1.3, 1.5, 1.8, 2.0, 3.0)),

                          # Input: Select the mu0 ----
                          selectInput(inputId = "p_mu0",
                                      label = "Mean:",
                                      choices = list_param_precomupteddataset$mu0),

                          # Input: Select the dispersion ----
                          selectInput(inputId = "p_dispersion",
                                      label = "Dispersion:",
                                      choices = list_param_precomupteddataset$dispersion),

                          # Input: Check box to choose which model to run ----
                          selectInput(inputId = "p_modelcheckGroup",
                                      label = "Models to run",
                                      choices = list("CytoGLMM - GLMM" = "cytoglmm",
                                                     "diffcyt - limma with fixed effect" = "testDS_limma_fixed",
                                                     "diffcyt - limma with random effect" = "testDS_limma_random",
                                                     "diffcyt - lme4" = "testDS_lmm")),

                          # Input: Action button to submit the parameters
                          actionButton(inputId = "p_goButton",
                                       label = "Get power"),

                          textOutput("na_error")


                        ),
                        # Main panel
                        mainPanel(
                          # textOutput("na_error"),
                          conditionalPanel("input.p_nb_donors == 'NA'",
                                           plotOutput("p_plot_pwr_nbdonors")
                          ),
                          conditionalPanel("input.p_nb_cells == 'NA'",
                                           plotOutput("p_plot_pwr_nbcells")
                          ),
                          conditionalPanel("input.p_subeff == 'NA'",
                                           plotOutput("p_plot_pwr_subject_effect")
                          ),
                          conditionalPanel("input.p_total_nb_markers == 'NA'",
                                           plotOutput("p_plot_pwr_total_nb_markers")
                          ),
                          conditionalPanel("input.p_mu0 == 'NA'",
                                           plotOutput("p_plot_pwr_mu0")
                          ),
                          conditionalPanel("input.p_dispersion == 'NA'",
                                           plotOutput("p_plot_pwr_dispersion")
                          ),
                          conditionalPanel("input.p_nb_donors != 'NA' &&
                                           input.p_nb_cells != 'NA' &&
                                           input.p_subeff != 'NA' &&
                                           input.p_total_nb_markers != 'NA' &&
                                           input.p_mu0 != 'NA' &&
                                           input.p_dispersion != 'NA'",
                                           textOutput(outputId = "p_power")
                          ))
                        # textOutput(outputId = "p_power"))
                      )
             ),
             # Personalized datasets
             tabPanel("Personalized dataset",
                      "This panel allows to compute the power for a chosen set of
                      parameter. The data are generated on request and it takes some time
                      to get the results.",
                      # Sidebar panel
                      sidebarLayout(

                        # Sidebar panel for inputs ----
                        sidebarPanel(

                          # Input: Field to enter the number of donors ----
                          numericInput(inputId = "nb_donors",
                                       label = "Number of donors:",
                                       value = 10,
                                       min = 3),
                          textOutput("nb_donors_warning"),

                          # Input: Slider for the subject effect ----
                          sliderInput(inputId = "subeff",
                                      label = "Subject effect:",
                                      min = 0.01,
                                      max = 3,
                                      value = 0.5),

                          # Input: Field to enter the number of cells per sample ----
                          numericInput(inputId = "nb_cells",
                                       label = "Number of cells per sample:",
                                       value = 500,
                                       min = 10),

                          # Input: Slider for the number of markers ----
                          sliderInput(inputId = "nb_tot_markers",
                                      label = "Total number of markers:",
                                      min = 3,
                                      max = 50,
                                      value = 20),

                          uiOutput("matrix"),

                          # Input: Check box to choose which model to run ----
                          checkboxGroupInput(inputId = "modelcheckGroup",
                                             label = "Models to run",
                                             choices = list("CytoGLMM - GLMM" = "cytoglmm",
                                                            # "CytoGLMM - GLM with bootstrap" = "cytoglm",
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
                        mainPanel(

                          # Output: Tabset w/ plot, summary, and table ----
                          tabsetPanel(type = "tabs",
                                      tabPanel("Summary",
                                               h4(textOutput(outputId = "title_selected_var")),
                                               uiOutput(outputId = "selected_var"),
                                               textOutput(outputId = "DE_names")),

                                      tabPanel("Power",
                                               plotOutput(outputId = "plot_pwr_eff"),
                                               h4(textOutput(outputId = "title_power")),
                                               DT::dataTableOutput(outputId = "results_pwr_with_effsize")),

                                      tabPanel("Details",
                                               # Observed variance and effect size
                                               h4(textOutput(outputId = "title_var")),
                                               DT::dataTableOutput(outputId = "tab_obs_variance"),
                                               ## Results
                                               # plotOutput(outputId = "plot_pval_dist"),
                                               DT::dataTableOutput(outputId = "results_models"))
                          )
                        )

                      )
             )
  )

)
