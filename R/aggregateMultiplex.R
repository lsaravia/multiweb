#' Aggregate a Multiplex Network into a Single Layer
#'
#' Combines a list of `igraph` layers into one aggregate network.
#' Edges present in multiple layers will have their weights summed.
#'
#' @param g_list A list of `igraph` objects (directed or undirected).
#' @param directed Logical; whether to treat the result as directed (default: `TRUE`).
#'
#' @return An `igraph` object representing the aggregate network.
#'
#' @import igraph
#' @export
aggregate_multiplex_network <- function(g_list, directed = TRUE) {
  if (!all(sapply(g_list, igraph::is.igraph))) stop("All elements of g_list must be igraph objects.")

  # Combine all node names
  all_nodes <- unique(unlist(lapply(g_list, function(g) V(g)$name)))

  # Create empty graph with all nodes
  g_aggr <- make_empty_graph(directed = directed) %>%
    add_vertices(length(all_nodes), name = all_nodes)

  # Merge edges across layers
  edge_df_list <- lapply(g_list, function(g) {
    el <- igraph::as_data_frame(g, what = "edges")
    if (is.null(el$weight)) el$weight <- 1
    el
  })

  combined_edges <- do.call(rbind, edge_df_list)

  # Sum weights of duplicate edges
  edge_sum <- combined_edges %>%
    group_by(from, to) %>%
    summarise(weight = sum(weight), .groups = "drop")

  # Add edges to the aggregate graph
  g_aggr <- add_edges(g_aggr, t(as.matrix(edge_sum[, c("from", "to")])))
  E(g_aggr)$weight <- edge_sum$weight

  return(g_aggr)
}
