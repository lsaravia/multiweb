% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/models.r
\name{generate_niche}
\alias{generate_niche}
\title{Generate Niche Model Food Web(s)}
\usage{
generate_niche(S, C, nsim = 1)
}
\arguments{
\item{S}{Integer. The number of species in the community. Must be \code{S > 1}.}

\item{C}{Numeric. The connectance (fraction of realized links). Must be \verb{0 < C ≤ 1}.}

\item{nsim}{Integer. The number of networks to generate (\verb{nsim ≥ 1}).}
}
\value{
If \code{nsim = 1}, returns a \strong{binary adjacency matrix (\verb{S × S})}.
If \code{nsim > 1}, returns a \strong{list of \code{igraph} objects}, each representing a food web.
}
\description{
This function generates one or multiple food webs using the \strong{Niche Model} proposed by Williams & Martinez (2000).
The model assumes that each species has a niche value and consumes resources within a defined range.
}
\examples{
generate_niche(20, 0.1)        # Single adjacency matrix
generate_niche(20, 0.1, nsim=5) # List of 5 food webs as igraph objects
}
\references{
Williams, R. J., and N. D. Martinez. 2000. Simple rules yield complex food webs. \emph{Nature} 404:180–183.
}
