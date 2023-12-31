% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_models_onesim.R
\name{function_to_compute_model_computation_onesimulation_modelchoice}
\alias{function_to_compute_model_computation_onesimulation_modelchoice}
\title{Run DS tests for one simulation.}
\usage{
function_to_compute_model_computation_onesimulation_modelchoice(
  onevariation,
  model
)
}
\arguments{
\item{onevariation}{list, of simulated data (output of the
function_wrapper_apply_simulation_nbtimes function).}

\item{model}{vector, name(s) of models to test.}
}
\value{
data.frame of results for each simulation all models combined.
}
\description{
Run DS tests for one simulation.
}
\details{
Function to run the DS tests through the different models.
}
\keyword{internal}
