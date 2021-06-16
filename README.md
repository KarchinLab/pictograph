# PICTograph

PICTograph is a computational tool developed to infer the clonal evolution of tumors from multi-region sequencing data. It models uncertainty in assigning mutations to subclones using a Bayesian hierarchical model and reduces the space of possible evolutionary trees by using constraints based on principles of sample-presence, lineage precedence, and sum condition. Highly probable evolutionary relationships that are recovered in multiple models can be highlighted using ensemble-based visualizations.

## Installation
```
devtools::install_github("KarchinLab/pictograph")
```
