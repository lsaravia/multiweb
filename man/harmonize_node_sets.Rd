% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_multiplex.R
\name{harmonize_node_sets}
\alias{harmonize_node_sets}
\title{Harmonize Node Sets Across Layers in a Multiplex Network}
\usage{
harmonize_node_sets(g.list)
}
\arguments{
\item{g.list}{A list of igraph objects representing the multiplex layers.}
}
\value{
A list of igraph objects with identical node sets.
}
\description{
This function ensures that all igraph layers in a multiplex network have the same set of nodes.
Nodes missing in any layer are added as isolated nodes (no edges).
}
\examples{
g.list <- harmonize_node_sets(g.list)
}
