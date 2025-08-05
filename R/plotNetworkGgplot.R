#' Plot an Ecological Network Organized by Trophic Level
#'
#' This function visualizes a food web using `ggplot2`, where nodes represent species and edges represent interactions.
#' Nodes are positioned by trophic level on the y-axis, and optionally, community modules on the x-axis.
#'
#' @param g An `igraph` object representing the ecological network.
#' @param vertexSizeFactor Numeric factor to scale node size (default: 5).
#' @param vertexSizeMin Minimum node size (default: 5).
#' @param modules Logical; if `TRUE`, nodes are grouped by community modules (default: `FALSE`).
#' @param weights Edge weights for community detection (default: `NA`).
#' @param node_weights Optional numeric vector for node weights. If `NULL` and the graph has a vertex attribute called `weight`, that attribute is used. If not, node degree is used as fallback.
#' @param community_obj Optional community detection object.
#' @param use_numbers Logical; if `TRUE`, label nodes by numeric ID (default: `FALSE`).
#' @param label_size Label font size (default: 4).
#' @param arrow_size Arrowhead size (default: 0.15).
#' @param shorten_factor Factor to shorten edges for arrows (default: 0.005).
#'
#' @return A `ggplot` object visualizing the trophic structure of the network.
#'
#' @export
plot_troph_level_ggplot <- function(
    g,
    vertexSizeFactor = 5,
    vertexSizeMin = 5,
    modules = FALSE,
    weights = NA,
    node_weights = NA,
    community_obj = NULL,
    use_numbers = FALSE,
    label_size = 4,
    arrow_size = 0.15,
    shorten_factor = 0.005
) {
  # --- Determine node weights ---
  # Node weights
  if (any(is.na(node_weights))) {
    node_weights <- degree(g, mode = "all")
  } else if (is.null(node_weights)) {
    if(!is.null(V(g)$weight)) {
       node_weights <- V(g)$weight
    } else {
      stop("`node_weights` must be provided or `V(g)$weights` must exist.")
    }
  }

  adj <- as_adjacency_matrix(g, sparse = FALSE)
  tl <- NetIndices::TrophInd(adj)

  if (is.null(V(g)$name)) {
    V(g)$name <- as.character(1:vcount(g))
  }

  nodes <- tibble(
    id = V(g)$name,
    weight = node_weights,
    size = log10(node_weights) * vertexSizeFactor + vertexSizeMin,
    trophic_level = tl$TL
  )

  if (use_numbers) {
    nodes <- nodes %>% mutate(numeric_id = row_number())
    label_column <- "numeric_id"
    node_map <- nodes %>% select(numeric_id, id)
  } else {
    label_column <- "id"
    node_map <- NULL
  }

  # X positions
  if (modules) {
    if (!is.null(community_obj)) {
      m <- community_obj
    } else {
      if (count_components(g) > 1) {
        if (!is.named(g)) V(g)$name <- as.character(1:vcount(g))
        dg <- components(g)
        V(g)$membership <- 0
        for (comp in unique(dg$membership)) {
          g1 <- induced_subgraph(g, which(dg$membership == comp))
          m <- cluster_spinglass(g1, weights = weights)
          if (length(m$membership) == 0) m$membership <- 1
          V(g)[V(g1)$name]$membership <- m$membership + max(V(g)$membership)
        }
        m$membership <- V(g)$membership
      } else {
        m <- cluster_spinglass(g, weights = weights)
      }
    }
    nodes <- nodes %>% mutate(module = factor(m$membership))
    nodes <- nodes %>% mutate(x = jitter(as.numeric(module), amount = 0.5))
  } else {
    nodes <- nodes %>% mutate(x = runif(n(), min = -1, max = 1))
  }

  nodes <- nodes %>% mutate(y = jitter(trophic_level, amount = 0.05))

  edges <- as.data.frame(as_edgelist(g)) %>%
    rename(from = V1, to = V2) %>%
    left_join(nodes, by = c("from" = "id")) %>%
    rename(x_from = x, y_from = y, size_from = size) %>%
    left_join(nodes, by = c("to" = "id")) %>%
    rename(x_to = x, y_to = y, size_to = size)

  module_lines <- NULL
  if (modules) {
    module_positions <- sort(unique(as.numeric(nodes$module)))
    if (length(module_positions) > 1) {
      line_positions <- head(module_positions, -1) + diff(module_positions) / 2
      module_lines <- tibble(x = line_positions)
    }
  }

  directed <- is_directed(g)
  if (directed) {
    edges <- edges %>%
      mutate(
        dx = x_to - x_from,
        dy = y_to - y_from,
        length = sqrt(dx^2 + dy^2),
        x_from = x_from + shorten_factor * dx / length * size_from,
        y_from = y_from + shorten_factor * dy / length * size_from,
        x_to = x_to - shorten_factor * dx / length * size_to,
        y_to = y_to - shorten_factor * dy / length * size_to
      )
  }

  arrow_style <- if (directed) arrow(type = "closed", length = unit(arrow_size, "inches")) else NULL

  p <- ggplot() +
    geom_segment(data = edges, aes(x = x_from, y = y_from, xend = x_to, yend = y_to),
                 color = "gray50", alpha = 0.6, size = 0.4, arrow = arrow_style) +
    geom_point(data = nodes, aes(x = x, y = y, size = size, color = trophic_level), alpha = 0.9) +
    geom_text_repel(data = nodes, aes(x = x, y = y, label = .data[[label_column]]),
                    max.overlaps = 15, size = label_size) +
    scale_color_viridis_c(option = "D") +
    theme_bw() +
    theme(
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(y = "Trophic Level")

  if (!is.null(module_lines)) {
    p <- p + geom_vline(data = module_lines, aes(xintercept = x), color = "black",
                        linetype = "dashed", alpha = 0.7)
  }

  if (modules) {
    p <- p + labs(x = "Community Module")
  } else {
    p <- p + theme(axis.title.x = element_blank())
  }

  return(p)
}


#' Plot an Ecological Network Organized by Trophic Level using ggraph
#'
#' This function visualizes a food web using `ggplot2`, where nodes represent species and edges represent interactions.
#' Nodes are positioned by trophic level on the y-axis, and optionally, community modules on the x-axis.
#' Is similar to [plot_troph_level_ggplot()], but uses `ggraph` for layout and rendering, which make better arrows.
#'
#' @param g An `igraph` network.
#' @param vertexSizeFactor Factor for node size scaling (default 5).
#' @param vertexSizeMin Minimum node size (default 5).
#' @param modules If TRUE, uses `cluster_spinglass` for community detection.
#' @param weights Edge weights for clustering (default NA).
#' @param community_obj Optional community detection object.
#' @param use_numbers Label nodes by numeric ID.
#' @param label_size Text size.
#' @param arrow_size Arrowhead size.
#' @param node_weights Optional numeric vector for node weights. If NA, uses degree.
#'        If it is NULL and the input graph has a ‘weight’ vertex attribute, then that attribute will be used.
#'
#' @return A `ggraph` plot.
#' @export
#'
#' @import igraph
#' @import ggraph
#' @import NetIndices
#' @import dplyr
#' @import ggrepel
#' @importFrom tibble tibble
#' @importFrom ggplot2 labs theme_bw theme element_blank element_text scale_color_viridis_c
plot_troph_level_ggraph <- function(g,
                                    vertexSizeFactor = 5,
                                    vertexSizeMin = 5,
                                    modules = FALSE,
                                    weights = NA,
                                    community_obj = NULL,
                                    use_numbers = FALSE,
                                    label_size = 4,
                                    arrow_size = 0.1,
                                    node_weights = NA) {
  # Compute trophic level
  adj <- as_adjacency_matrix(g, sparse = FALSE)
  tl <- NetIndices::TrophInd(adj)

  if (is.null(V(g)$name)) {
    V(g)$name <- as.character(1:vcount(g))
  }

  # Node weights
  if (any(is.na(node_weights))) {
    V(g)$weight <- degree(g, mode = "all")
  } else if (!is.null(node_weights)) {
    V(g)$weight <- node_weights
  } else if(is.null(V(g)$weight)){
    stop("`node_weights` must be provided or `V(g)$weights` must exist.")
  }

  # Node sizes
  V(g)$size<- V(g)$weight * vertexSizeFactor
  V(g)$size<- ifelse(V(g)$size < vertexSizeMin, vertexSizeMin, V(g)$size)
  V(g)$trophic_level <- tl$TL

  # Modules
  if (modules) {
    if (!is.null(community_obj)) {
      m <- community_obj
    } else {
      m <- cluster_spinglass(g, weights = weights)
    }
    V(g)$module <- factor(m$membership)
  } else {
    V(g)$module <- factor(1)  # dummy for layout
  }

  # Numeric labels if needed
  if (use_numbers) {
    V(g)$numeric_id <- seq_len(vcount(g))
  }

  # Layout coordinates: x = module, y = trophic level
  if (modules) {
    x_pos <- jitter(as.numeric(V(g)$module), amount = 0.5)
  } else {
    x_pos <- runif(vcount(g), -1, 1)
  }
  y_pos <- V(g)$trophic_level

  layout_df <- tibble::tibble(x = x_pos, y = y_pos)

  # Build plot
  p <- ggraph(g, layout = "manual", x = layout_df$x, y = layout_df$y) +
    geom_edge_link(
      aes(),
      arrow = arrow(type = "closed", length = unit(arrow_size, "inches")),
      end_cap = circle(arrow_size, 'inches'),
      start_cap = circle(arrow_size, 'inches'),
      edge_colour = "grey50",
      alpha = 0.6
    ) +
    geom_node_point(
      aes(size = size, color = trophic_level),
      alpha = 0.9
    ) +
    geom_text_repel(
      aes(x = layout_df$x, y = layout_df$y,label = if (use_numbers) numeric_id else name),
      size = label_size,
      max.overlaps = 20
    ) +
    scale_color_viridis_c() +
    theme_bw() +
    labs(y = "Trophic Level") +
    theme(
      legend.position = "none",
      panel.grid = element_blank()
    )

  # Add module lines
  if (modules) {
    module_positions <- sort(unique(as.numeric(V(g)$module)))
    if (length(module_positions) > 1) {
      line_positions <- head(module_positions, -1) + diff(module_positions) / 2
      p <- p + geom_vline(
        xintercept = line_positions,
        linetype = "dashed",
        color = "black",
        alpha = 0.7
      ) +
        labs(x = "Community Module")
    }
  } else {
    p <- p + theme(axis.title.x = element_blank())
  }

  return(p)
}

