% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/modularity.R
\name{run_infomap_multi}
\alias{run_infomap_multi}
\title{Run Infomap on Multilayer Networks with Intralayer Links}
\usage{
run_infomap_multi(
  igraph_list,
  layer_names = NULL,
  infomap_path = "infomap",
  output_dir = tempdir(),
  directed = TRUE,
  two_level = TRUE,
  multilayer_relax_rate = 0.15,
  seed = 123
)
}
\arguments{
\item{igraph_list}{A list of \code{igraph} objects, each representing a network layer.}

\item{layer_names}{A character vector with layer names corresponding to each igraph object, default NULL.}

\item{infomap_path}{Path to the Infomap binary (default assumes it's in the system \code{PATH}).}

\item{output_dir}{A temporary directory for Infomap results. Default is \code{tempdir()}.}

\item{directed}{Logical: Treat the network as directed (default \code{TRUE}).}

\item{two_level}{Logical: Use two-level Infomap (default \code{TRUE}).}

\item{multilayer_relax_rate}{Numeric: Relaxation rate for multilayer links (default \code{0.15}).}

\item{seed}{Numeric: Random seed for Infomap (default \code{123}).}
}
\value{
A data frame containing the detected modules with columns:
\item{module}{Module/community assignment from Infomap.}
\item{node}{The actual node name from the original igraph objects.}
\item{layer}{The corresponding layer name from \code{layer_names}.}
\item{flow}{The fraction of flow assigned to the module.}
}
\description{
This function exports a multilayer network (as a list of igraph objects) into Infomap's required format,
runs Infomap via an external binary, and imports the detected communities back into R as a data frame.
It requires Infomap to be installed on the system from \href{https://www.mapequation.org/infomap/#Install}{Infomap}.
}
\details{
If the network has the \code{weight} attribute, it will be used as the weight of the edges.
}
\examples{
\dontrun{
  ig_list <- list(ig_p1, ig_p2, ig_p3)  # List of igraph objects
  layer_names <- c("Pond 1", "Pond 2", "Pond 3")  # Corresponding layer names
  result <- run_infomap_multi(ig_list, layer_names)
  print(result)
}

}
