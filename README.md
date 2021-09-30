# CyTOFpower

## Description

This package is a tool to predict the power of CyTOF experiments, in the context 
of differential state analyses using `CytoGLMM` or `diffcyt` packages. 
This package provides a shiny app with two options to predict the power of an 
experiment. One part generates in-sicilico CyTOF data, using users input parameters 
(for instance the number of markers, the number of cells, or the expected fold 
change for the differential expressed markers). And the other part allows the user 
to browse in a grid of parameters for which the power was precomputed.

## Example

More information are available in the vignette:

```	r
vignette("CyTOFpower", package = "CyTOFpower")
```
