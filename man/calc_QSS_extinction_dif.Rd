% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calcExtinctionIndices.r
\name{calc_QSS_extinction_dif}
\alias{calc_QSS_extinction_dif}
\title{Calculates the QSS difference between the full network and the network minus
one species}
\usage{
calc_QSS_extinction_dif(g, sp_list, nsim = 1000, ncores = 4, istrength = FALSE)
}
\arguments{
\item{g}{igraph network}

\item{sp_list}{list with the species/nodes we will delete for the comparison}

\item{nsim}{number of simulations to calculate QSS}

\item{ncores}{number of cores used to perform the operation}

\item{istrength}{if TRUE takes the weight attribute of the network as interaction strength to
calculate QSS.}
}
\value{
a data.frame with nsim rows for each species deleted with:
\itemize{
\item the node deleted
\item The QSS of the complete network for nsim simulations
\item The QSS of the network with the deleted node for nsim simulations
\item The difference between the two previous QSS
}
}
\description{
The Quasi-sign stability is estimated with of the community matrix (Jacobian) and characterizes
the linear stability of the network. This uses the function \code{\link[=calc_QSS]{calc_QSS()}} so it can take
into account the interaction strength if weights are present.
}
\examples{
\dontrun{
g <- netData[[1]]

# Generate random weights
#
V(g)$weight <-  runif(vcount(g))

# Without interaction strength
#
calc_QSS_extinction_dif(g,V(g)$name[1:3],nsim=10,istrength = FALSE)

# With interaction strength
#
calc_QSS_extinction_dif(g,V(g)$name[1:3],nsim=10,istrength = TRUE)
}
}
