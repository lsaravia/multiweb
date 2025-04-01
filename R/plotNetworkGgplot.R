#' Plot an Ecological Network Organized by Trophic Level
#'
#' This function visualizes a food web using `ggplot2`, where nodes represent species and edges represent interactions.
#' Nodes are positioned by trophic level on the y-axis, and optionally, community modules on the x-axis.
#'
#' @param g An `igraph` object representing the ecological network.
#' @param vertexSizeFactor Numeric factor to determine node size based on degree (default: 5).
#' @param vertexSizeMin Numeric value for the minimum node size (default: 5).
#' @param modules Logical; if `TRUE`, nodes are grouped by community modules detected using `cluster_spinglass()` (default: `FALSE`).
#' @param weights A numeric vector, `NULL`, or `NA`, specifying edge weights for community detection.
#'   If `NULL`, the function uses the `weight` attribute of `g` if available. If `NA`, weights are ignored (default: `NA`).
#' @param community_obj An optional community detection object. If `NULL`, the function calculates modules automatically (default: `NULL`).
#' @param maxTL Numeric, maximum trophic level to display on the y-axis. If `NULL`, it is determined automatically (default: `NULL`).
#' @param use_numbers Logical; if `TRUE`, nodes are labeled with numeric IDs instead of species names to reduce label overlap in large networks (default: `FALSE`).
#' @param label_size Numeric, font size for node labels (default: 4).
#' @param arrow_size Numeric, size of the arrowheads for directed edges (default: 0.15).
#'
#' @return A `ggplot` object visualizing the trophic structure of the network.
#'   If `use_numbers = TRUE`, it also returns a `tibble` mapping numeric labels to species names.
#'
#' @importFrom igraph degree get.adjacency V count_components cluster_spinglass induced_subgraph components
#' @importFrom NetIndices TrophInd
#' @importFrom ggplot2 ggplot geom_segment geom_point theme_bw theme element_blank element_text labs
#' @importFrom ggrepel geom_text_repel
#' @importFrom dplyr mutate row_number select left_join
#' @importFrom tibble tibble as_tibble
#' @export
#'
#' @examples
#'
#' g <- netData$Angola
#' plot_troph_level_ggplot(g, modules = TRUE, use_numbers = TRUE)
#'
#'
plot_troph_level_ggplot <- function(g, vertexSizeFactor = 5, vertexSizeMin = 5, modules = FALSE,
                                    weights = NA, community_obj = NULL, maxTL = NULL,
                                    use_numbers = FALSE, label_size = 4, arrow_size = 0.15) {

  # Compute node degree and trophic level
  deg <- degree(g, mode = "all")
  adj <- as_adjacency_matrix(g, sparse = FALSE)
  tl <- TrophInd(adj)

  # Ensure node names exist, otherwise assign numeric IDs
  if (is.null(V(g)$name)) {
    V(g)$name <- as.character(1:vcount(g))  # Assign numbers as names
  }
  # Create node tibble
  nodes <- tibble(
    id = V(g)$name,
    degree = deg,
    size = log10(deg + 1) * vertexSizeFactor + vertexSizeMin,
    trophic_level = tl$TL
  )

  # Assign numbers to nodes if `use_numbers = TRUE`
  if (use_numbers) {
    nodes <- nodes %>% mutate(numeric_id = row_number())
    label_column <- "numeric_id"
    node_map <- nodes %>% select(numeric_id, id)  # Store mapping of numbers to names
  } else {
    label_column <- "id"
    node_map <- NULL
  }

  # Compute layout (x = modularity, y = trophic level)
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

  # Convert edges to a tibble
  edges <- as.data.frame(as_edgelist(g)) %>%
    rename(from = V1, to = V2) %>%
    left_join(nodes, by = c("from" = "id")) %>%
    rename(x_from = x, y_from = y, size_from = size) %>%
    left_join(nodes, by = c("to" = "id")) %>%
    rename(x_to = x, y_to = y, size_to = size)


  # Create vertical separator lines if modules are enabled
  module_lines <- NULL
  if (modules) {
    module_positions <- sort(unique(as.numeric(nodes$module)))
    if (length(module_positions) > 1) {
      line_positions <- head(module_positions, -1) + diff(module_positions) / 2  # Midpoints
      module_lines <- tibble(x = line_positions)
    }
  }

  # Check if the graph is directed
  directed <- is_directed(g)

  if (directed) {
    # Adjust arrow placement by shortening the edges
    shorten_factor <- 0.005  # Adjust this factor if needed
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


  # Define arrow style with adjustable arrow size
  arrow_style <- if (directed) arrow(type = "closed", length = unit(arrow_size, "inches")) else NULL

  # Create base ggplot network plot
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


  # Add vertical separator lines **only between modules**
  if (!is.null(module_lines)) {
    p <- p + geom_vline(data = module_lines, aes(xintercept = x), color = "black",
                        linetype = "dashed", alpha = 0.7)
  }

  # Add x-axis label **only if modules are used**
  if (modules) {
    p <- p + labs(x = "Community Module")
  } else {
    p <- p + theme(axis.title.x = element_blank())  # Remove x-axis label
  }

  print(p)

  # Return node mapping if using numbers
  if (use_numbers) {
    return(node_map)
  }
}
