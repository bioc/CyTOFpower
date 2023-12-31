# Data for shiny interface
precomputed_data <- readRDS(system.file("extdata",
                                        "df_precomputed_datasets.rds",
                                        package = "CyTOFpower"))
list_param_precomupteddataset <-list(
  "nb_donor" = c(sort(unique(precomputed_data$nb_donor)), NA),
  "total_nb_markers" = c(sort(unique(precomputed_data$total_nb_markers)), NA),
  "nb_DEmarkers" = sort(unique(precomputed_data$nb_DEmarkers)),
  "subject_effect" = c(sort(unique(precomputed_data$subject_effect)), NA),
  "mu0" = c(sort(unique(precomputed_data$mu0)), NA),
  "dispersion" = c(sort(unique(precomputed_data$dispersion)), NA),
  "nb_cell_per_sample" = c(sort(unique(precomputed_data$nb_cell_per_sample)), NA)
)

#' Shiny app
#'
#' Interactive shiny app to predict the power of a CyTOF experiment.
#'
#' @return Interactive shiny app.
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
