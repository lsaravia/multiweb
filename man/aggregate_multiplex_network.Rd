% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aggregateMultiplex.R
\name{aggregate_multiplex_network}
\alias{aggregate_multiplex_network}
\title{Aggregate a Multiplex Network into a Single Layer}
\usage{
aggregate_multiplex_network(g_list, directed = TRUE)
}
\arguments{
\item{g_list}{A list of \code{igraph} objects (directed or undirected).}

\item{directed}{Logical; whether to treat the result as directed (default: \code{TRUE}).}
}
\value{
An \code{igraph} object representing the aggregate network.
}
\description{
Combines a list of \code{igraph} layers into one aggregate network.
Edges present in multiple layers will have their weights summed.
}
