% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_models.R
\name{function_formula_contrast_diffcytDSLMM_randomeffect}
\alias{function_formula_contrast_diffcytDSLMM_randomeffect}
\title{Formula and contrast matrix for diffcyt-DS-LMM with random effect.}
\usage{
function_formula_contrast_diffcytDSLMM_randomeffect(df_experiment_info)
}
\arguments{
\item{df_experiment_info}{data.frame, information about the experiment (donor IDs,
group IDs and sample IDs).}
}
\value{
list with 3 slots:
    - formula: the formula for the model;
    - contrast: the contrast matrix.
}
\description{
Formula and contrast matrix for diffcyt-DS-LMM with random effect.
}
\details{
Function to create formula and contrast for diffcyt-DS-LMM with random effect.
}
\keyword{internal}
