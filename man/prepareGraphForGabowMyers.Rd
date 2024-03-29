% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gabowmyers.R
\name{prepareGraphForGabowMyers}
\alias{prepareGraphForGabowMyers}
\title{Create tibble of possible edges from CCF values based on sample-presence}
\usage{
prepareGraphForGabowMyers(w, chains, input_data)
}
\arguments{
\item{w}{matrix of CCF values (rows = clusters, columns = samples)}

\item{chains}{list of MCMC chains (must contain w_chain and z_chain)}

\item{input_data}{list of input data; same as supplied for clustering}
}
\value{
graph_G tibble of possible edges with columns edge, parent, child
}
\description{
Create tibble of possible edges from CCF values based on sample-presence
}
