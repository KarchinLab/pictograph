% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/subclone-proportions.R
\name{forceCCFs}
\alias{forceCCFs}
\title{Force CCFs to comply with lineage precedence and sum condition}
\usage{
forceCCFs(w_mat, tree_edges)
}
\arguments{
\item{w_mat}{Matrix of CCF estimates (from \code{estimateCCFs})}

\item{tree_edges}{Tibble of tree edges with columns edge, parent, and child}
}
\description{
Force CCFs to comply with lineage precedence and sum condition
}
