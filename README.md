# CyTOFpower

## Description

This package is a tool to predict the power of CyTOF 
experiments, where differential state analyses are performed, using `CytoGLMM` or
`diffcyt` packages. This package provides a shiny app with two of checking for 
the power. One part generates in-sicilico CyTOF data, using users input parameters 
(for instance the number of markers, the number of cells, or the expected fold 
change for the differential expressed markers). And the other part allows the user 
to browse in a grid of parameters where the power was precomputed.

## Example

More information are available in the vignette:

```	r
vignette("CyTOFpower", package = "CyTOFpower")
```
