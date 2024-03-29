% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calcTopologicalIndices.r
\name{calc_modularity}
\alias{calc_modularity}
\title{Calculation of Modularity for a list of igraph objects}
\usage{
calc_modularity(ig, ncores = 0, weights = NA)
}
\arguments{
\item{ig}{list of igraph objects to calculate modularity}

\item{ncores}{number of cores used to compute in parallel, if 0 sequential processing is used.}

\item{weights}{The weights of the edges. Either a numeric vector or NULL or NA. If it is null and the input graph has a ‘weight’ edge attribute
then that will be used. If NULL and no such attribute is present then the edges will have equal weights.
Set this to NA if the graph was a ‘weight’ edge attribute, but you don't want to use it for community detection.}
}
\value{
a data.frame with the field Modularity
}
\description{
The function calculates modularity of the list of networks in the nullDist parameter.
Modularity is calculated using the \code{\link[igraph:cluster_spinglass]{igraph::cluster_spinglass()}}
Only works for one component networks.
}
\examples{
\dontrun{
nullg <- generateERbasal(netData[[1]],10)
calc_modularity(nullg)
}

}
