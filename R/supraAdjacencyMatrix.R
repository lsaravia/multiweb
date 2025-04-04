#' Generate Node aligned Interlayer Edges for a Multilayer Network from the function `convert_to_intra_format`
#'
#' @param intra_data A list with `$vertices` (node info) and `$intra` (intralayer edges).
#' @param interlayer_weight Numeric, the weight for interlayer links (default: 1).
#' @return A dataframe with interlayer edges.
add_interlayer_links <- function(intra_data, interlayer_weight = 1) {
  nodes <- intra_data$vertices
  intra_edges <- intra_data$intra

  # Extract unique node-layer combinations from both columns
  node_layers <- unique(rbind(
    data.frame(layer_id = intra_edges$layer_id, node_id = intra_edges$node_id1),
    data.frame(layer_id = intra_edges$layer_id, node_id = intra_edges$node_id2)
  ))

  # Identify nodes that exist in multiple layers
  multilayer_nodes <- table(node_layers$node_id)
  aligned_nodes <- as.integer(names(multilayer_nodes[multilayer_nodes > 1]))

  # Generate interlayer edges only for nodes in multiple layers
  inter_edges <- do.call(rbind, lapply(aligned_nodes, function(node) {
    layers <- sort(unique(node_layers$layer_id[node_layers$node_id == node]))

    if (length(layers) > 1) {
      df <- data.frame(
        node_id = rep(node, length(layers) - 1),
        layer_from = head(layers, -1),
        layer_to = tail(layers, -1),
        weight = interlayer_weight
      )

      # Create bidirectional edges by swapping from ↔ to
      df_bidirectional <- df[, c("node_id", "layer_to", "layer_from", "weight")]
      colnames(df_bidirectional) <- colnames(df)

      return(rbind(df, df_bidirectional))  # Combine both directions
    }
    return(NULL)
  }))

  return(inter_edges)
}



#' Convert Multilayer Network to Supra-Adjacency Matrix (SAM)
#'
#' @param igraph_list A list of `igraph` objects (each representing a layer).
#' @param layer_names (Optional) A vector of names for each layer.
#' @param interlayer_weight Numeric, the weight for interlayer links (default: 1).
#' @param isDirected Logical, whether the network is directed.
#' @param sparse Logical, whether to return the supra-adjacency matrix as a sparse matrix (default: TRUE).
#' @param use_names Logical; if `TRUE`, edges will use node names instead of numeric IDs (default: `FALSE`).
#' @param interlayer Logical; if `FALSE`, no interlayer links are added (default: `TRUE`).
#' @param clean Logical; if `TRUE`, removes empty rows/columns (state nodes with no edges) from the output (default: TRUE).
#' @return A list with:
#'   \item{supra_matrix}{Supra-adjacency matrix with optional inter-layer links.}
#'   \item{state_nodes_map}{Mapping of state nodes (layer-node combinations).}
#' @import igraph
#' @import Matrix
#' @export
convert_to_supra_adjacency <- function(igraph_list, layer_names = NULL, interlayer_weight = 1,
                                       isDirected = TRUE, sparse = TRUE, use_names = FALSE,
                                       interlayer = TRUE, clean = TRUE) {
  # Step 1: Convert to intra-layer format
  intra_data <- convert_to_intra_format(igraph_list)
  num_layers <- length(igraph_list)
  num_nodes <- nrow(intra_data$vertices)

  # Assign default layer names if not provided
  if (is.null(layer_names)) {
    layer_names <- paste0("layer_", seq_len(num_layers))
  }

  # Step 2: Create state-node mapping
  state_nodes_map <- expand.grid(
    node_id = intra_data$vertices$node_id,
    layer_id = seq_len(num_layers)
  )
  state_nodes_map <- merge(state_nodes_map, intra_data$vertices, by = "node_id", all.x = TRUE)
  state_nodes_map <- merge(state_nodes_map, data.frame(layer_id = seq_len(num_layers), layer_name = layer_names),
                           by = "layer_id", all.x = TRUE)

  # Assign unique state-node IDs
  state_nodes_map <- state_nodes_map[order(state_nodes_map$layer_id, state_nodes_map$node_id), ]
  state_nodes_map$sn_id <- seq_len(nrow(state_nodes_map))

  # Use node names if requested
  if (use_names) {
    state_nodes_map$tuple <- paste(state_nodes_map$layer_name, state_nodes_map$name, sep = "_")
  } else {
    state_nodes_map$tuple <- paste(state_nodes_map$layer_name, state_nodes_map$node_id, sep = "_")
  }

  # Step 3: Generate inter-layer edges (if enabled)
  interlayer_edges <- if (interlayer) add_interlayer_links(intra_data, interlayer_weight) else NULL

  # Step 4: Convert intra-layer edges to extended edge list format
  intra_edges <- intra_data$intra
  extended_intra <- data.frame(
    node_id1 = intra_edges$node_id1,
    layer_id1 = intra_edges$layer_id,
    node_id2 = intra_edges$node_id2,
    layer_id2 = intra_edges$layer_id,
    weight = intra_edges$weight
  )

  # Step 5: Convert inter-layer edges (if available)
  extended_inter <- if (interlayer) {
    data.frame(
      node_id1 = interlayer_edges$node_id,
      layer_id1 = interlayer_edges$layer_from,
      node_id2 = interlayer_edges$node_id,
      layer_id2 = interlayer_edges$layer_to,
      weight = interlayer_edges$weight
    )
  } else NULL

  # Merge intra-layer and inter-layer edges
  full_edge_list <- if (!is.null(extended_inter)) rbind(extended_intra, extended_inter) else extended_intra

  # Step 6: Convert to Supra-Adjacency Matrix (SAM)
  SAM <- BuildSupraAdjacencyMatrixFromExtendedEdgelist(
    mEdges = full_edge_list,
    Layers = num_layers,
    Nodes = num_nodes,
    isDirected = isDirected
  )

  # Assign row/column names if using node names
  if (use_names) {
    rownames(SAM) <- state_nodes_map$tuple
    colnames(SAM) <- state_nodes_map$tuple
  }

  # Step 7: Remove empty rows/columns if requested
  if (clean) {
    active <- Matrix::rowSums(SAM) + Matrix::colSums(SAM) > 0
    SAM <- SAM[active, active, drop = FALSE]
    state_nodes_map <- state_nodes_map[active, , drop = FALSE]
  }

  return(list(supra_matrix = SAM, state_nodes_map = state_nodes_map))
}

#' Build supra-adjacency matrix from edge lists
#'
#'
#' @param mEdges data frame with extended edge list, i.e a data frame with
#'   columns: \code{node.from, layer.from, node.to, layer.to weight}
#' @param Layers scalar, number of layers
#' @param Nodes scalar, number of nodes
#' @param isDirected logical
#' @return
#' Supra-adjacency matrix, a square matrix of dimension \code{Nodes * Layers}.
BuildSupraAdjacencyMatrixFromExtendedEdgelist <-
  function(mEdges, Layers, Nodes, isDirected) {

  if (max(max(mEdges[, 2]), max(mEdges[, 4])) != Layers) {
    stop("Error: expected number of layers does not match the data. Aborting process.")
  }

  edges <- data.frame(
    from = mEdges[, 1] + Nodes * (mEdges[, 2] - 1),
    to = mEdges[, 3] + Nodes * (mEdges[, 4] - 1),
    weight = mEdges[, 5]
  )

  M <-
    Matrix::sparseMatrix(
      i = edges$from,
      j = edges$to,
      x = edges$weight,
      dims = c(Nodes * Layers, Nodes * Layers)
    )

  if (sum(abs(M - Matrix::t(M))) > 1e-12 && isDirected == FALSE) {
    message(
      "WARNING: The input data is directed but isDirected=FALSE, I am symmetrizing by average."
    )
    M <- (M + Matrix::t(M)) / 2
  }
  return(M)
}
