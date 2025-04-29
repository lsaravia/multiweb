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
#'
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
  edges <- as.data.frame(as_edgelist(graph))
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
  codelength_line <- clu_lines[grep("^# codelength", clu_lines)]
  codelength <- if (length(codelength_line) > 0) {
    as.numeric(sub(".*codelength ([0-9\\.]+) bits.*", "\\1", codelength_line))
  } else {
    NA
  }

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

#' Convert a List of igraph Objects to Intralayer Edge Format
#'
#' This function converts a multilayer network, represented as a list of igraph objects,
#' into a standardized format suitable for Infomap and other multilayer network analyses.
#'
#' Each node is assigned a unique numeric ID, and edges across layers are recorded
#' with layer-specific identifiers.
#'
#' @param igraph_list A list of `igraph` objects, each representing a network layer.
#' @param use_names Logical; if `TRUE`, edges will use node names instead of numeric IDs (default: `FALSE`).
#' @return A list with two data frames:
#'   \item{vertices}{A data frame containing node_id and corresponding node names.}
#'   \item{intra}{A data frame containing intra-layer edges with columns: `layer_id`, `node_id1`, `node_id2`, `weight`.}
#'
#' @examples
#' \dontrun{
#'   multi_format <- convert_to_intra_format(list(layer1, layer2))
#'   print(multi_format$vertices)
#'   print(multi_format$intra)
#' }
#'
#' @import igraph
#' @export
convert_to_intra_format <- function(igraph_list, use_names = FALSE) {
  # Extract unique node names from all layers
  all_nodes <- unique(unlist(lapply(igraph_list, function(g) V(g)$name)))
  node_df <- data.frame(node_id = seq_along(all_nodes), name = all_nodes)

  # Create a mapping from node name to numeric ID
  node_map <- setNames(node_df$node_id, node_df$name)

  # Collect edges across layers
  edge_list <- do.call(rbind, lapply(seq_along(igraph_list), function(i) {
    g <- igraph_list[[i]]
    layer_id <- i
    edges <- igraph::as_data_frame(g, what = "edges")

    # Ensure weight column exists
    if (!"weight" %in% colnames(edges)) {
      edges$weight <- 1  # Default weight if none exists
    }

    # Map node names to numeric IDs
    tibble(
      layer_id = layer_id,
      node_id1 = if (use_names) edges$from else node_map[edges$from],
      node_id2 = if (use_names) edges$to else node_map[edges$to],
      weight = edges$weight
    )
  }))

  return(list(vertices = node_df, intra = edge_list))
}

#' Write a Multilayer Network to a File in Infomap-Compatible Format
#'
#' This function writes a multilayer network (represented as a list of igraph objects)
#' to a text file in Infomap's required *Intra format.
#'
#' The output file consists of:
#' - A list of unique vertices with their numeric IDs.
#' - Intra-layer edges with weights.
#'
#' @param igraph_list A list of `igraph` objects, each representing a layer.
#' @param file_path A string specifying the path to the output file.
#' @return A list with two data frames:
#'   \item{vertices}{A data frame of node IDs and names.}
#'   \item{intra}{A data frame of intra-layer edges with columns: `layer_id`, `node_id1`, `node_id2`, `weight`.}
#'
#' @examples
#' \dontrun{
#'   write_multilayer_network(list(layer1, layer2), "network.net")
#' }
#'
#' @import igraph
#' @export
write_multilayer_network <- function(igraph_list, file_path) {
  # Convert the igraph objects to a multilayer dataframe
  multilayer_df <- convert_to_intra_format(igraph_list)
  # Extract vertices and intra-layer edges correctly
  vertices <- multilayer_df$vertices
  intra_edges <- multilayer_df$intra

  # Open the file for writing
  con <- file(file_path, "w")

  # Write header
  writeLines("# A multilayer network using *Intra format", con)

  # Write vertices
  writeLines(paste0("*Vertices ", nrow(vertices)), con)
  writeLines("# node_id name", con)
  write.table(vertices, con, row.names = FALSE, col.names = FALSE, quote = TRUE)

  # Write intra-layer edges
  writeLines("*Intra", con)
  writeLines("# layer_id node_id node_id weight", con)
  write.table(intra_edges, con, row.names = FALSE, col.names = FALSE, quote = FALSE)

  # Close the file
  close(con)
  return(multilayer_df)
}

#' Run Infomap on Multilayer Networks with Intralayer Links
#'
#' This function exports a multilayer network (as a list of igraph objects) into Infomap's required format,
#' runs Infomap via an external binary, and imports the detected communities back into R as a data frame.
#' It requires Infomap to be installed on the system. See installation instructions at
#' [Infomap](https://www.mapequation.org/infomap/#Install).
#'
#' If the network has the `weight` attribute, it will be used as the weight of the edges.
#' Since the exported multilayer network structure assumes only `*Intra` links (i.e., no explicit inter-layer links),
#' Infomap internally constructs inter-layer connections for each physical node using the `--multilayer-relax-rate`
#' parameter (default: 0.15). This creates inter-layer links between all instances of the same node across layers.
#' The weights of these links are proportional to the weighted out-degree of the node in the target layer,
#' resulting in non-uniform inter-layer transition probabilities.
#'
#' This implicit construction of inter-layer links enables Infomap to model node-aligned multiplex dynamics
#' without requiring an explicit `*Inter` section in the input file.

#'
#' @param igraph_list A list of `igraph` objects, each representing a network layer.
#' @param layer_names A character vector with layer names corresponding to each igraph object (default: `NULL`).
#' @param infomap_path Character string specifying the path to the Infomap binary (default: `"infomap"`, assumes it is in the system `PATH`).
#' @param output_dir Character string specifying the directory for Infomap results (default: `tempdir()`).
#' @param directed Logical; if `TRUE`, the network is treated as directed (default: `TRUE`).
#' @param two_level Logical; if `TRUE`, uses two-level Infomap instead of hierarchical Infomap (default: `TRUE`).
#' @param multilayer_relax_rate Numeric; relaxation rate for multilayer links (default: `0.15`).
#' @param seed Numeric; random seed for Infomap to ensure reproducibility (default: `123`).
#'
#' @return A data frame containing the detected modules with columns:
#' \describe{
#'   \item{module}{Module/community assignment from Infomap.}
#'   \item{node}{The node name from the original igraph objects.}
#'   \item{layer}{The corresponding layer name from `layer_names`.}
#'   \item{flow}{The fraction of flow assigned to the module.}
#' }
#'
#' @references
#' D. Edler, A. Holmgren and M. Rosvall, *The MapEquation software package*, available online at \url{https://www.mapequation.org}.
#'
#' @importFrom igraph as_edgelist E V
#' @importFrom dplyr mutate select rename left_join bind_rows
#'
#' @examples
#' \dontrun{
#' # Load network data
#' fileName <- system.file("extdata", package = "multiweb")
#' dn <- list.files(fileName, pattern = "^Kefi2015.*\\.txt$")
#' g <- readNetwork(dn, fileName, skipColumn = 2)
#' class(g)
#' names(g) <- c("Negative", "Positive", "Trophic")
#'
#' # Run Infomap
#' run_infomap_multi(g)
#' }
#'
#' @export
run_infomap_multi <- function(igraph_list, layer_names=NULL, infomap_path = "infomap", output_dir = tempdir(),
                              directed = TRUE, two_level = TRUE, multilayer_relax_rate = 0.15, seed = 123) {

  # Define file paths
  net_file <- file.path(output_dir, "network.net")
  clu_file <- file.path(output_dir, "network.clu")
  clu_states_file <- file.path(output_dir, "network_states.clu")

  # Export graph to Infomap's required format
  ml_df <- write_multilayer_network(igraph_list, net_file)

  # Construct Infomap command
  options <- "--silent --clu"
  if (directed) options <- paste(options, "-d")
  if (two_level) options <- paste(options, "-2")
  options <- paste(options, "--multilayer-relax-rate", multilayer_relax_rate)
  options <- paste(options, "--seed", seed)

  # Run Infomap
  command <- sprintf('%s "%s" "%s" %s', infomap_path, net_file, output_dir, options)
  system(command, intern = TRUE)

  # Read Infomap output
  clu_lines <- readLines(clu_file)

  # Extract codelength from output
  codelength_line <- clu_lines[grep("^# codelength", clu_lines)]
  codelength <- if (length(codelength_line) > 0) {
    as.numeric(sub(".*codelength ([0-9\\.]+) bits.*", "\\1", codelength_line))
  } else {
    NA
  }

  # Read the states.clu file
  clu_data <- read.delim(clu_states_file, comment.char = "#", header = FALSE, sep = " ",
                         col.names = c("state_id", "module", "flow", "node_id", "layer_id"))

  # Ensure correct order
  clu_data <- clu_data[order(clu_data$layer_id, clu_data$node_id), ]

  # Use ml_df$vertices to replace numeric node_id with actual node names
  clu_data <- clu_data %>%
    left_join(ml_df$vertices, by = c("node_id")) %>%
    rename(node = name) %>%
    dplyr::select(-state_id, -node_id)

  # Replace layer_id with layer_name if names are provided
  if (!is.null(layer_names) && length(layer_names) >= max(clu_data$layer_id, na.rm = TRUE)) {
    clu_data$layer <- layer_names[clu_data$layer_id]
  } else {
    clu_data <- clu_data %>% rename(layer = layer_id)  # Keep numeric ID if names are missing
  }


  result <- list(
    communities = clu_data %>% dplyr::select(module, node, layer, flow),
    codelength = codelength
  )
  return(result)
}


