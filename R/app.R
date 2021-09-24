# Data for shiny interface
precomputed_data <- utils::read.table(system.file("extdata",
                                        "df_precomputed_datasets.txt",
                                        package = "CyTOFpower"))
list_param_precomupteddataset <-list(
  "nb_donor" = c(unique(precomputed_data$nb_donor), NA),
  "total_nb_markers" = c(unique(precomputed_data$total_nb_markers), NA),
  "nb_DEmarkers" = unique(precomputed_data$nb_DEmarkers),
  "subject_effect" = c(unique(precomputed_data$subject_effect), NA),
  "mu0" = c(unique(precomputed_data$mu0), NA),
  "dispersion" = c(unique(precomputed_data$dispersion), NA),
  "nb_cell_per_sample" = c(unique(precomputed_data$nb_cell_per_sample), NA)
)

#' Shiny app
#'
#' Interactive shiny app to predict the power of a CyTOF experiment.
#'
#' @return Shiny app.
#' @examples
#' # Launch the shiny app
#' if (interactive()) {
#'   CyTOFpower()
#' }
#' @export
CyTOFpower <- function(){
  app <- shinyApp(ui = appUI, server = appServer)
  runApp(app)
}
