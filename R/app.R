# Data for shiny interface
precomputed_data <- utils::read.table(system.file("extdata",
                                        "df_precomputed_datasets.txt",
                                        package = "CyTOFpower"))
list_param_precomupteddataset <-list(
  "nb_donor" = c(unique(precomputed_data$nb_donor), NA), #c(5, 10, 20, 30, NA), #c(3, 5, 10, 20, 30, 40, 50, NA),
  "total_nb_markers" = c(unique(precomputed_data$total_nb_markers), NA), #c(10, 20, 30, 40, NA), #c(10, 20, 30, 40, 50, 60),
  "nb_DEmarkers" = unique(precomputed_data$nb_DEmarkers), #7,
  "subject_effect" = c(unique(precomputed_data$subject_effect), NA), #c(0.1, 0.5, NA),  #c(seq(0.1, 1, by = 0.2), NA),
  "mu0" = c(unique(precomputed_data$mu0), NA), #c(seq(1, 10, by = 2), NA), #c(1:10, NA),
  "dispersion" = c(unique(precomputed_data$dispersion), NA), #c(seq(0.5, 4, by = 1), NA), #by = 0.5
  "nb_cell_per_sample" = c(unique(precomputed_data$nb_cell_per_sample), NA)  #c(1000, 5000, 10000, 20000, NA)  #c(500, 1000, 2000, 5000, 10000, 15000, 20000, NA)
)

#' Shiny app
#'
#' Interactive shiny app to predict the power of a CyTOF experiment.
#'
#' @export
CyTOFpower <- function(){
  app <- shinyApp(ui = appUI, server = appServer)
  runApp(app)
}
