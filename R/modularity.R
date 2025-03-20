#' Calculation of Modularity for a List of igraph Objects
#'
#' This function calculates the modularity of a list of networks provided in the `ig` parameter.
#' Modularity measures the strength of division of a network into modules (communities).
#'
#' **Note:** `cluster_spinglass()` only works on networks with a single connected component.
#'
#' @param ig A list of igraph objects for which modularity will be calculated.
#' @param ncores Number of cores used for parallel computation. If set to 0, sequential processing is used.
#' @param weights Edge weights. Can be a numeric vector, `NULL`, or `NA`:
#'   - If `NULL`, the function checks for a 'weight' edge attribute and uses it if present.
#'   - If `NA`, weights are ignored even if the graph has a 'weight' edge attribute.
#' @param cluster_function The **igraph** clustering function to use, or a function that returns modularity.
#'   default `cluster_spinglass`.
#'
#' @return A data frame with a column `Modularity` containing the modularity values for each network.
#'
#' @export
#'
#' @import igraph
#' @importFrom future.apply future_lapply
#' @importFrom future sequential multisession
#' @importFrom dplyr bind_rows
#'
#' @examples
#' \dontrun{
#' nullg <- generateERbasal(netData[[1]], 10)
#' calc_modularity(nullg, cluster_function = cluster_infomap)
#' }
#'
calc_modularity <- function(ig, ncores = 0, cluster_function = cluster_spinglass) {

  # Ensure input is a list of igraph objects
  if (inherits(ig, "igraph")) {
    ig <- list(ig)
  } else if (!inherits(ig[[1]], "igraph")) {
    stop("The parameter 'ig' must be a list of igraph objects.")
  }

  # Register parallel backend
  if (ncores > 0) {
    cn <- future::availableCores()
    if (ncores > cn) ncores <- cn
    future::plan(multisession, workers = ncores)
    on.exit(future::plan(sequential))
  } else {
    future::plan(sequential)
  }

  # Function to compute modularity (handling different weight parameter names) CHECK
  compute_modularity <- function(g) {

    community <- cluster_function(g)
    data.frame(Modularity = modularity(community))
  }

  # Compute modularity in parallel
  df <- future_lapply(ig, compute_modularity, future.seed = TRUE)

  return(bind_rows(df))
}

#' Interface R and igraph with Infomap (external binary)
#'
#' This function exports an igraph network to Infomap's required format, runs Infomap,
#' and imports the detected communities back into R as an igraph-compatible cluster object.
#' It requires the installation of infomap in the system from https://www.mapequation.org/infomap/#Install
#' if the network has the weight attribute, it will be used as the weight of the edges.
#'
#' @param graph An igraph object.
#' @param infomap_path Path to the Infomap binary (default assumes it's in system PATH).
#' @param output_dir Temporary directory for Infomap results.
#' @param directed Boolean: Treat the network as directed (default TRUE).
#' @param two_level Boolean: Use two-level Infomap (default TRUE).
#' @param seed Numeric: Random seed for Infomap (default 123).
#' @return An igraph community object similar to igraph's cluster_infomap.
#'
#' @examples
#' # Example with the intermediate files in actual folder
#' py_infomap <- run_infomap(netData[[1]],output_dir=".", weighted = FALSE)
#' results <- calculate_eigencentrality(A)
#' membership(py_infomap)
#' modularity(py_infomap)
#'
#' @import igraph
#' @export
run_infomap <- function(graph, infomap_path = "infomap", output_dir = tempdir(),
                        directed = TRUE, two_level = TRUE, seed = 123) {

  # Define file paths
  net_file <- file.path(output_dir, "network.net")
  clu_file <- file.path(output_dir, "network.clu")

  # Export graph to Infomap's own format
  edges <- as.data.frame(get.edgelist(graph))
  colnames(edges) <- c("source", "target")

  if(is.null(edge_attr(graph, "weight"))) {
    edges$weight <- 1  # Default weight if unweighted
  } else {
    edges$weight <- E(graph)$weight
  }

  # Ensure node IDs are consecutive integers (Infomap requires this)
  node_ids <- setNames(seq_along(V(graph)), V(graph)$name)
  edges$source <- node_ids[edges$source]
  edges$target <- node_ids[edges$target]

  # Write Infomap-compatible file
  writeLines(c("# Directed network with weights",
               paste(edges$source, edges$target, edges$weight)), net_file)

  # Construct Infomap command
  options <- "--silent --clu"
  if (directed) options <- paste(options, "-d")
  if (two_level) options <- paste(options, "-2")
  options <- paste(options, "--seed", seed)

  # Run Infomap
  command <- sprintf('%s "%s" "%s" %s', infomap_path, net_file, output_dir, options)
  system(command, intern = TRUE)

  # Read Infomap output
  clu_lines <- readLines(clu_file)

  # Extract codelength from output
  codelength_line <- clu_lines[grep("codelength", clu_lines)]
  if (length(codelength_line) > 0) {
    codelength_match <- regmatches(codelength_line, regexpr("[0-9]+\\.[0-9]+", codelength_line))
    codelength <- as.numeric(codelength_match)
  } else {
    codelength <- NA  # Assign NA if no match is found
  }
  #codelength <- as.numeric(sub(".*codelength ([0-9.]+) bits.*", "\\1", codelength_line))[1]

  # Read the .clu file

  clu_data <- read.delim(clu_file, comment.char = "#", header = FALSE, sep = " ", col.names = c("node_id", "module", "flow"))
  # Ensure correct order by sorting by node_id
  clu_data <- clu_data[order(clu_data$node_id), ]

  # Assign communities to igraph object
  membership <- clu_data$module

  community_obj <- make_clusters(graph, membership = membership, algorithm = "infomap")
  community_obj$codelength <- codelength

  # Assign names if available
  if (!is.null(V(graph)$name)) {
    community_obj$names <- V(graph)$name
  }
  return(community_obj)
}

