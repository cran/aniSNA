% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/regression_analysis.R
\name{regression_slope_analyze}
\alias{regression_slope_analyze}
\title{To perform regression analysis for local network metrics}
\usage{
regression_slope_analyze(
  network,
  n_simulations = 10,
  subsampling_proportion = c(0.1, 0.3, 0.5, 0.7, 0.9),
  network_metrics = c("degree", "strength", "betweenness", "clustering_coefficient",
    "eigenvector_centrality")
)
}
\arguments{
\item{network}{An igraph graph object consisting of observed network}

\item{n_simulations}{Number of sub-samples to be obtained at each level}

\item{subsampling_proportion}{A vector depicting proportions of sub-sampled nodes}

\item{network_metrics}{A vector depicting names of local global network metrics}
}
\value{
A list of network metrics of class list_regression_matrices. Each element of list is a matrix whose columns 
        correspond to subsampling_proportion and rows correspond to n_simulations.
        The entries of the matrix provide value of the slope of regression when the 
        nodal values in sub-sampled network are regressed upon the values of the same 
        nodes in the full network for the corresponding metric.
}
\description{
To perform regression analysis for local network metrics
}
\examples{
\donttest{
data(elk_network_2010)
regression_slope_analyze(elk_network_2010)
}
}
