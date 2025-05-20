#' 3D Visualization of a Multiplex Ecological Network
#'
#' This function generates a 3D visualization of a multiplex network composed of multiple layers
#' (e.g., trophic, mutualistic, and competitive networks) using `rgl`. Layers are displayed stacked
#' along the Z-axis. Node positions are shared across layers, and nodes can be colored by module
#' membership or functional type. An optional aggregate layer (merging all layers) can also be displayed.
#'
#' @param g.list A list of `igraph` objects, each representing a layer of the multiplex network.
#' @param layer.colors A vector of colors corresponding to each network layer.
#' @param as.undirected Logical; if `TRUE`, converts each layer to undirected for visualization purposes (default: `TRUE`).
#' @param layer.layout A layout matrix shared by all layers, or `"auto"` to use `layoutMultiplex()` (default: `"auto"`).
#' @param layer.labels A character vector for labeling layers. `"auto"` assigns "Layer 1", "Layer 2", etc. (default: `"auto"`).
#' @param layer.labels.cex Numeric; text size for layer labels (default: `2`).
#' @param edge.colors Edge colors (a single color or vector per layer), or `"auto"` to match layer color (default: `"auto"`).
#' @param edge.normalize Logical; if `TRUE`, normalizes edge weights using logarithmic scale (default: `FALSE`).
#' @param edge.size.scale Numeric or vector to control edge thickness (default: `1`).
#' @param node.colors `"auto"` to use layer colors, or an `N x L` matrix for node colors per layer.
#' @param node.size.values `"auto"` to size nodes by strength, or a fixed numeric value (default: `0.5`).
#' @param node.size.scale Numeric or vector to scale node size (default: `1`).
#' @param node.alpha Transparency for nodes (0–1; default: `1`).
#' @param edge.alpha Transparency for edges (0–1; default: `1`).
#' @param layer.alpha A numeric vector with layer transparency, or `"auto"` to set 0.5 for all layers.
#' @param layout A character string indicating the layout method for `layoutMultiplex()` (default: `"fr"`).
#' @param show.nodeLabels Logical; if `TRUE`, node labels are shown (default: `FALSE`).
#' @param show.aggregate Logical; if `TRUE`, a layer with the aggregated network is shown (default: `FALSE`).
#' @param aggr.alpha Transparency for the aggregate layer (default: `"auto"` which is set to 0.5).
#' @param aggr.color Color of the background layer for the aggregate network (default: `"#dadada"`).
#' @param node.colors.aggr Node color used in the aggregate network (default: `"#dadada"`).
#' @param layer.scale Scale factor for the extent of each layer (default: `2`).
#' @param layer.shift.x Horizontal X shift per layer (default: `0`).
#' @param layer.shift.y Horizontal Y shift per layer (default: `0`).
#' @param layer.space Vertical spacing between layers along Z (default: `1.5`).
#' @param FOV Field of view for the `rgl` perspective visualization (default: `30`).
#'
#' @details
#' Each layer is visualized as a semi-transparent plane, on which the respective network is plotted.
#' Nodes are shared across layers and maintain consistent positions. Optionally, the aggregated network
#' (union of all layers) can be visualized at the top. This function is designed for ecological multiplex
#' networks and facilitates interpretation of modular structure, interaction types, or perturbation effects.
#'
#' This function is adapted from the `plot_multiplex3D()` function
#' in the [`muxViz`](https://manlius.github.io/muxViz/) project.
#'
#' @import igraph rgl
#' @export
#'
#' @examples
#' \dontrun{
#' library(multiweb)
#' fileName <- system.file("extdata", package = "multiweb")
#' dn <- list.files(fileName, pattern = "^Kefi2015.*\\.txt$")
#' g_list <- readNetwork(dn, fileName, skipColumn = 2)
#' names(g_list) <- c("Negative", "Positive", "Trophic")
#'
#'
#' # Visualize with community-colored nodes
#' plot_multi3D(
#'   g.list = g_list,
#'   layer.colors = RColorBrewer::brewer.pal(3, "Set2"),
#'   node.colors = "auto",
#'   node.size.scale = 0.5,
#'   node.size.values = "auto",
#'   node.colors.aggr = "#999999",
#'   show.aggregate = TRUE
#' )
#' }
plot_multi3D <- function (g.list, layer.colors, as.undirected = T, layer.layout = "auto",
          layer.labels = "auto", layer.labels.cex = 2, edge.colors = "auto",
          edge.normalize = F, edge.size.scale = 1, node.colors = "auto",
          node.size.values = 0.5, node.size.scale = 1, node.alpha = 1,
          edge.alpha = 1, layer.alpha = "auto", layout = "fr", show.nodeLabels = F,
          show.aggregate = F, aggr.alpha = "auto", aggr.color = "#dadada",
          node.colors.aggr = "#dadada", layer.scale = 2, layer.shift.x = 0,
          layer.shift.y = 0, layer.space = 1.5, FOV = 30)
{
  mypal <- layer.colors
  Layers <- length(g.list)
  Nodes <- igraph::vcount(g.list[[1]])
  if (!is.matrix(layer.layout) && layer.layout == "auto") {
    lay <- layoutMultiplex(g.list, layout = layout, ggplot.format = F,
                           box = T)
  }
  else {
    lay <- layer.layout
  }
  if (identical(layer.alpha,"auto")) {
    layer.alpha <- rep(0.5, Layers)
  }
  if (is.na(layer.labels[1]) || is.null(layer.labels[1])) {
    layer.labels <- NA
  }
  else {
    if (layer.labels[1] == "auto" || length(layer.labels) !=
        Layers) {
      layer.labels <- paste("Layer", 1:Layers)
    }
    if (show.aggregate && (!is.na(layer.labels[1]) && !is.null(layer.labels[1]))) {
      layer.labels <- c(layer.labels, "Aggregate")
    }
  }
  if (length(node.size.scale) == 1) {
    node.size.scale <- rep(node.size.scale, Layers)
  }
  if (length(edge.size.scale) == 1) {
    edge.size.scale <- rep(edge.size.scale, Layers)
  }
  LAYER_SCALE <- layer.scale
  LAYER_SHIFT_X <- layer.shift.x
  LAYER_SHIFT_Y <- layer.shift.y
  LAYER_SPACE <- layer.space
  PLOT_FOV <- FOV
  d <- 0
  clear3d()
  bg3d(col = "white")
  for (l in 1:Layers) {
    if (as.undirected) {
      g.list[[l]] <- igraph::as_undirected(g.list[[l]])
    }
    if (identical(node.size.values,"auto")) {
      igraph::V(g.list[[l]])$size <- 3 * node.size.scale[l] *
        sqrt(igraph::strength(g.list[[l]]))
    }
    else {
      igraph::V(g.list[[l]])$size <- node.size.values[,l] *
        node.size.scale[l]
    }
    if (!is.matrix(node.colors)) {
      if (node.colors[1] == "auto") {
        node.col <- layer.colors[l]
      }
      else {
        node.col <- node.colors
      }
      igraph::V(g.list[[l]])$color <- node.col
    }
    else {
      igraph::V(g.list[[l]])$color <- node.colors[, l]
    }
    if (show.nodeLabels) {
      degs <- igraph::degree(g.list[[l]])
      labels <- seq_len(igraph::gorder(g.list[[l]]))  # Use numeric IDs
      igraph::V(g.list[[l]])$label <- ifelse(degs > 0, labels, NA)
    } else {
      igraph::V(g.list[[l]])$label <- NA
    }

    if (edge.colors == "auto") {
      edge.col <- layer.colors[l]
    }
    else {
      edge.col <- edge.colors
    }
    igraph::E(g.list[[l]])$color <- edge.col
    if (!is.null(igraph::E(g.list[[l]])$weight)) {
      igraph::E(g.list[[l]])$width <- igraph::E(g.list[[l]])$weight
    }
    else {
      igraph::E(g.list[[l]])$width <- 1
    }
    if (edge.normalize) {
      igraph::E(g.list[[l]])$width <- edge.size.scale[l] *
        log(1 + igraph::E(g.list[[l]])$width)/max(log(1 +
                                                        igraph::E(g.list[[l]])$width))
    }
    if (show.aggregate) {
      d <- -1 + LAYER_SCALE * LAYER_SPACE * l/(Layers +
                                                 1)
    }
    else {
      d <- -1 + LAYER_SCALE * LAYER_SPACE * l/Layers
    }
    layout.layer <- matrix(0, nrow = Nodes, ncol = 3)
    layout.layer[, 1] <- lay[, 1] + (l - 1) * LAYER_SHIFT_X
    layout.layer[, 2] <- lay[, 2] + (l - 1) * LAYER_SHIFT_Y
    layout.layer[, 3] <- d
    x <- c(-1, -1, -1 + LAYER_SCALE, -1 + LAYER_SCALE) +
      (l - 1) * LAYER_SHIFT_X
    y <- c(-1 + LAYER_SCALE, -1, -1, -1 + LAYER_SCALE) +
      (l - 1) * LAYER_SHIFT_Y
    z <- c(d, d, d, d)
    quads3d(x, y, z, alpha = layer.alpha[[l]], col = layer.colors[[l]],
            add = T)
    igraph::rglplot(g.list[[l]], layout = layout.layer,
                    rescale = F)
    if (!is.na(layer.labels[1]) && !is.null(layer.labels[1])) {
      text3d(-1 + (l - 1) * LAYER_SHIFT_X, -1 + (l - 1) *
               LAYER_SHIFT_Y, d + 0.1, text = layer.labels[l],
             adj = 0.2, color = "black", family = "sans",
             cex = layer.labels.cex)
    }
  }
  if (show.aggregate) {
    g.aggr <- GetAggregateNetworkFromNetworkList(g.list)
    if (identical(node.size.values,"auto")) {
      igraph::V(g.aggr)$size <- 3 * node.size.scale[l] *
        sqrt(igraph::strength(g.aggr))
    }
    else {
      igraph::V(g.aggr)$size <- sumR(node.size.values,2) * node.size.scale[l]
    }
    igraph::V(g.aggr)$color <- node.colors.aggr
    if (show.nodeLabels) {
      igraph::V(g.aggr)$label <- 1:igraph::gorder(g.aggr)
    }
    else {
      igraph::V(g.aggr)$label <- NA
    }
    igraph::E(g.aggr)$color <- aggr.color
    if (!is.null(igraph::E(g.aggr)$weight)) {
      igraph::E(g.aggr)$width <- igraph::E(g.aggr)$weight
    }
    else {
      igraph::E(g.aggr)$width <- 1
    }
    l <- Layers + 1
    d <- -1 + LAYER_SCALE * LAYER_SPACE * l/(Layers + 1)
    layout.layer <- matrix(0, nrow = Nodes, ncol = 3)
    layout.layer[, 1] <- lay[, 1] + (l - 1) * LAYER_SHIFT_X
    layout.layer[, 2] <- lay[, 2] + (l - 1) * LAYER_SHIFT_Y
    layout.layer[, 3] <- d
    x <- c(-1, -1, -1 + LAYER_SCALE, -1 + LAYER_SCALE) +
      (l - 1) * LAYER_SHIFT_X
    y <- c(-1 + LAYER_SCALE, -1, -1, -1 + LAYER_SCALE) +
      (l - 1) * LAYER_SHIFT_Y
    z <- c(d, d, d, d)
    if (aggr.alpha == "auto") {
      quads3d(x, y, z, alpha = 0.5, col = aggr.color,
              add = T)
    }
    else {
      quads3d(x, y, z, alpha = aggr.alpha, col = aggr.color,
              add = T)
    }
    igraph::rglplot(g.aggr, layout = layout.layer, rescale = F)
    if (!is.na(layer.labels[1]) && !is.null(layer.labels[1])) {
      text3d(-1 + (l - 1) * LAYER_SHIFT_X, -1 + (l - 1) *
               LAYER_SHIFT_Y, d + 0.1, text = "Aggregate",
             adj = 0.2, color = "black", family = "sans",
             cex = layer.labels.cex)
    }
  }
  M <- matrix(0, ncol = 4, nrow = 4)
  M[1, ] <- c(0.54, 0, 0.84, 0)
  M[2, ] <- c(0.33, 0.92, -0.22, 0)
  M[3, ] <- c(-0.77, 0.39, 0.5, 0)
  M[4, ] <- c(0, 0, 0, 1)
  par3d(FOV = PLOT_FOV, userMatrix = M)
}


#' Convert Infomap Community Output to muxViz-Compatible Format
#'
#' This function converts the community detection result from `run_infomap_multi()` (from the `multiweb` package)
#' into a format compatible with `muxViz`, allowing direct use with visualization function `plot_multimodules()`.
#'
#' It creates two key data frames:
#' \itemize{
#'   \item \code{membership.multi}: Community membership for each node-layer combination.
#'   \item \code{membership.aggr}: Aggregate community membership per physical node (across layers),
#'   based on the highest flow.
#' }
#'
#' @param communities A data frame returned by `run_infomap_multi()$communities`, containing columns `node`, `layer`, `module`, and `flow`.
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{membership.multi}{A data frame with columns `node`, `layer`, and `module` representing the state node assignment. Node and layer are returned as integers (1-based index).}
#'   \item{membership.aggr}{A data frame with columns `node` and `module`, where each physical node is assigned to its dominant module across all layers.}
#'   \item{modules.multi}{Integer; the number of unique modules detected in the multilayer network.}
#'   \item{modules.aggr}{Integer; the number of unique modules in the aggregated network.}
#' }
#'
#' @examples
#' \dontrun{
#' # Run Infomap community detection on a multiplex ecological network
#' g_list <- readNetwork(dn, fileName, skipColumn = 2)
#' names(g_list) <- c("Negative", "Positive", "Trophic")
#' result <- run_infomap_multi(g_list, layer_names = names(g_list))
#'
#' # Convert for muxViz plotting
#' mux_data <- convert_infomap_result_to_muxviz(result$communities)
#'
#' # Visualize
#' plot_multimodules(mux_data)
#' }
#'
#' @seealso \code{\link{run_infomap_multi}}, \code{\link{plot_multi3D}}, \code{\link{plot_multimodules}}
#' @export
convert_infomap_result_to_muxviz <- function(communities) {
  # Aseguramos el orden de capas y nodos
  layers <- unique(communities$layer)
  nodes <- unique(communities$node)

  # Crear membership.multi con índices enteros
  node_index <- setNames(seq_along(nodes), nodes)
  layer_index <- setNames(seq_along(layers), layers)

  membership.multi <- communities %>%
    mutate(
      node = node_index[node],
      layer = layer_index[layer]
    ) %>%
    dplyr::select(node, layer, module) %>%
    arrange(layer, node)

  # Obtener módulo dominante por nodo
  membership.aggr <- communities %>%
    group_by(node, module) %>%
    summarise(flow = sum(flow), .groups = "drop") %>%
    group_by(node) %>%
    slice_max(flow, n = 1) %>%
    ungroup() %>%
    mutate(node = node_index[node]) %>%
    dplyr::select(node, module)

  # Devolver en el formato esperado por muxViz
  return(list(
    membership.multi = membership.multi,
    membership.aggr = membership.aggr,
    modules.multi = length(unique(membership.multi$module)),
    modules.aggr = length(unique(membership.aggr$module))
  ))
}


#' Harmonize Node Sets Across Layers in a Multiplex Network
#'
#' This function ensures that all igraph layers in a multiplex network have the same set of nodes.
#' Nodes missing in any layer are added as isolated nodes (no edges).
#'
#' @param g.list A list of igraph objects representing the multiplex layers.
#' @return A list of igraph objects with identical node sets.
#'
#' @examples
#' g.list <- harmonize_node_sets(g.list)
#' @export
harmonize_node_sets <- function(g.list) {
  # Get union of all node names across all layers
  all_nodes <- sort(unique(unlist(lapply(g.list, function(g) V(g)$name))))

  # Add missing nodes to each layer
  g.list <- lapply(g.list, function(g) {
    current_nodes <- V(g)$name
    missing_nodes <- setdiff(all_nodes, current_nodes)

    if (length(missing_nodes) > 0) {
      g <- igraph::add_vertices(g, nv = length(missing_nodes), name = missing_nodes)
    }

    return(g)
  })

  return(g.list)
}

#' Plot Modular Structure of a Multiplex Ecological Network
#'
#' This function visualizes each layer of a multiplex ecological network, with nodes colored
#' by module membership and sized by Infomap flow values. It returns a list of `ggplot` objects.
#'
#' @param g_list A named list of `igraph` objects (e.g., from `readNetwork()`), each representing one layer.
#' @param communities A data frame from `run_infomap_multi()$communities`, containing `node`, `module`, `layer`, `flow`.
#' @param module_palette Optional named color palette for module IDs (default uses `RColorBrewer::Set2`).
#' @param layout_graph Optional igraph object for layout (default is aggregate union of all layers).
#' @param scale_factor Numeric multiplier for node size scaling by flow (default: 500).
#'
#' @return A named list of `ggplot` objects, one for each network layer.
#'
#' @examples
#' fileName <- system.file("extdata", package = "multiweb")
#' dn <- list.files(fileName, pattern = "^Kefi2015.*\\.txt$")
#' g_list <- readNetwork(dn, fileName, skipColumn = 2)
#' names(g_list) <- c("Negative", "Positive", "Trophic")
#'
#' res <- run_infomap_multi(g_list, layer_names = names(g_list))
#' plots <- plot_multiplex_modules(g_list, res$communities)
#' cowplot::plot_grid(plotlist = plots, ncol = 3)
#'
#' @import ggplot2 dplyr ggraph tidygraph igraph
#' @export
plot_multiplex_modules <- function(g_list,
                                   communities,
                                   module_palette = NULL,
                                   layout_graph = NULL,
                                   scale_factor = 500) {

  if (is.null(names(g_list))) stop("g_list must be a named list of igraph objects.")

  # Compute shared layout using union of all layers
  if (is.null(layout_graph)) {
    layout_graph <- GetAggregateNetworkFromNetworkList(g_list)
  }
  lay <- layout_with_fr(layout_graph)
  layout_df <- data.frame(name = V(layout_graph)$name, x = lay[, 1], y = lay[, 2])

  # Get module set and build palette
  modules <- sort(unique(communities$module))
  if (is.null(module_palette)) {
    module_palette <- setNames(RColorBrewer::brewer.pal(min(length(modules), 8), "Set2"), modules)
  }

  # Generate plots per layer
  plot_list <- lapply(names(g_list), function(lname) {
    g_layer <- g_list[[lname]]
    node_names <- V(g_layer)$name

    # Get module + flow info
    df_comm <- communities %>%
      filter(layer == lname) %>%
      group_by(node) %>%
      slice_max(order_by = flow, n = 1, with_ties = FALSE) %>%
      ungroup()

    # Match layout
    layout_layer <- layout_df %>% filter(name %in% node_names)

    # Build tidygraph
    g_tbl <- as_tbl_graph(g_layer) %>%
      left_join(df_comm, by = c("name" = "node")) %>%
      mutate(
        module = factor(module),
        size = ifelse(!is.na(flow), flow * scale_factor, 2)
      )

    # Return ggplot
    ggraph(g_tbl, layout = "manual", x = layout_layer$x, y = layout_layer$y) +
      geom_edge_link(color = "gray80", alpha = 0.4) +
      geom_node_point(aes(color = module, size = size)) +
      scale_color_manual(values = module_palette, na.translate = FALSE) +
      scale_size_continuous(range = c(1, 8)) +
      labs(title = lname) +
      theme_void() +
      theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        legend.position = "none"
      )
  })

  names(plot_list) <- names(g_list)
  return(plot_list)
}

