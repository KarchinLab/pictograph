# PICTograph

PICTograph is a computational tool developed to infer the clonal evolution of tumors from multi-region sequencing data. It models uncertainty in assigning mutations to subclones using a Bayesian hierarchical model and reduces the space of possible evolutionary trees by using constraints based on principles of sample-presence, lineage precedence, and sum condition. Highly probable evolutionary relationships that are recovered in multiple models can be highlighted using ensemble-based visualizations.

## Installation

PICTograph uses the JAGS library for Bayesian data analysis, which is installed outside of R. JAGS can be downloaded and installed for your OS [here](https://mcmc-jags.sourceforge.io/).

To install PICTograph from GitHub, start R and enter: 
```
devtools::install_github("KarchinLab/pictograph", build_vignettes = TRUE)
```

## Tutorial

Demo code on a toy example can be found in the vignette: 
```
library(pictograph)
vignette("pictograph", package = "pictograph")
```
