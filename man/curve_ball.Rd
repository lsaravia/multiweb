% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/models.r
\name{curve_ball}
\alias{curve_ball}
\alias{curveBall}
\title{Curbe ball algorithm to generate random networks with given an igraph network}
\usage{
curve_ball(g, nsim = 1000, istrength = FALSE)
}
\arguments{
\item{g}{igraph object to extract adjacency matrix}

\item{nsim}{number of generated random networks}

\item{istrength}{if TRUE the edge attribute weight is taken as interaction strength and randomized keeping the column sum equal to the original matrix.}
}
\value{
a list of randomized igraph objects
}
\description{
This extract the adjacency matrix from an igraph object and generates a randomized version
of the network with the same row and column totals, the results have the same in and out degree sequence
than the original network. The diagonals are not avoided so it can generate self-links or cannibalism in
the context of food-webs. In the case that the network has multiple components (or disconnected networks)
the algorithm simulates around the same number of components, if the original network has 1 component
the algorithm enforces that the results have all 1 component. If the edge attribute weight is present
and istrength=TRUE then the weigth is additionally randomized keeping the column sum equal to the original matrix.
}
\details{
Based on:
}
\examples{

curve_ball(netData[[1]])

}
\references{
Strona, G. et al. 2014. A fast and unbiased procedure to randomize ecological binary matrices with
fixed row and column totals. -Nat. Comm. 5: 4114. doi: 10.1038/ncomms5114
}
