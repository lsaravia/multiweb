#' Plot an Ecological Network Organized by Trophic Level using visNetwork
#'
#' This function visualizes a food web using `visNetwork`, where nodes represent species and edges represent interactions.
#' Nodes are positioned by trophic level on the y-axis, and optionally, community modules on the x-axis.
#' Unlike [plot_troph_level_ggplot()] and [plot_troph_level_ggraph()], this produces an interactive HTML widget
#' that supports dragging, zooming, hovering for details, and node search.
#'
#' @param g An `igraph` object representing the ecological network.
#' @param vertexSizeMin Minimum node size (default: 8).
#' @param vertexSizeMax Maximum node size (default: 25).
#' @param modules Logical; if `TRUE`, nodes are grouped by community modules along the x-axis (default: `FALSE`).
#' @param weights Edge weights used for community detection when `modules = TRUE` and `community_obj` is not
#'   supplied (default: `NA`).
#' @param community_obj Optional pre-computed community detection object (e.g. from `cluster_spinglass()`).
#'   If `NULL` and `modules = TRUE`, community detection is run internally via `cluster_spinglass()`.
#' @param node_weights Optional numeric vector used to size nodes. If `NA` (default), node degree is used.
#' @param module_spacing Horizontal spacing multiplier between community modules when `modules = TRUE`
#'   (default: 350). When `modules = FALSE`, node x-coordinates are instead spread across a wide random
#'   range scaled by `module_spacing * 1.5`.
#' @param physics Layout mode, one of `"x"` or `"full"` (default: `"x"`).
#'   \itemize{
#'     \item `"x"`: nodes can move horizontally but their vertical position stays fixed at their
#'       trophic level; physics remains active continuously to allow ongoing horizontal reflow.
#'     \item `"full"`: nodes move freely in both x and y under a force-directed layout
#'       (`forceAtlas2Based`). A small `centralGravity` anchors the system to prevent unbounded drift,
#'       and physics is automatically disabled once stabilization finishes to avoid continuous
#'       rotation/jitter of the final layout.
#'   }
#' @param highlight Logical; if `TRUE`, hovering over a node highlights its nearest neighbors (default: `TRUE`).
#' @param search Logical; if `TRUE`, adds a node-search dropdown to the widget (default: `TRUE`).
#' @param label_size Font size, in pixels, used for node labels (default: 16). Apparent label size can
#'   differ between `physics` modes even at the same value, since `"full"` typically stabilizes into a
#'   more spread-out layout that is zoomed further out to fit the view; adjust this parameter if labels
#'   look too large or too small for a given mode.
#'
#' @return A `visNetwork` HTML widget visualizing the trophic structure of the network.
#'
#' @export
#'
#' @importFrom igraph as_adjacency_matrix V degree cluster_spinglass membership
#'   vertex_attr_names vcount ends E
#' @importFrom NetIndices TrophInd
#' @importFrom visNetwork visNetwork visNodes visEdges visOptions visInteraction
#'   visPhysics visEvents
#' @importFrom viridisLite viridis
plot_troph_level_visNetwork <- function(
    g,
    vertexSizeMin = 8,
    vertexSizeMax = 25,
    modules = FALSE,
    weights = NA,
    community_obj = NULL,
    node_weights = NA,
    module_spacing = 350,
    physics = c("x","full"),
    highlight = TRUE,
    search = TRUE,
    label_size = 16
){

  physics <- match.arg(physics)

  ##------------------------------------------------------------
  ## Trophic level
  ##------------------------------------------------------------

  adj <- as_adjacency_matrix(g, sparse = FALSE)
  tl <- NetIndices::TrophInd(adj)

  V(g)$trophic_level <- tl$TL

  ##------------------------------------------------------------
  ## Node weights
  ##------------------------------------------------------------

  if(length(node_weights)==1 && is.na(node_weights)){

    node_weights <- degree(g, mode="all")

  }

  V(g)$weight <- node_weights

  ##------------------------------------------------------------
  ## Community modules
  ##------------------------------------------------------------

  if(modules){

    if(is.null(community_obj))
      community_obj <- cluster_spinglass(g, weights = weights)

    V(g)$module <- membership(community_obj)

  }else{

    V(g)$module <- 1

  }

  ##------------------------------------------------------------
  ## Labels
  ##------------------------------------------------------------

  labels <-
    if("label" %in% vertex_attr_names(g))
      V(g)$label
  else
    V(g)$name

  ##------------------------------------------------------------
  ## Coordinates
  ##------------------------------------------------------------
  ## Si hay módulos, se agrupan las x por módulo (con jitter).
  ## Si no hay módulos, se distribuyen las x en un rango amplio
  ## para que los nodos no queden amontonados en una franja angosta.

  x <- if(modules){

    jitter(
      as.numeric(V(g)$module),
      amount = .45
    )

  }else{

    runif(vcount(g), min = -1, max = 1)

  }

  y <- -V(g)$trophic_level

  ##------------------------------------------------------------
  ## Node size
  ##------------------------------------------------------------

  s <- sqrt(V(g)$weight)

  if(diff(range(s)) > 0){

    s <- (s-min(s))/(max(s)-min(s))

    s <- vertexSizeMin +
      s*(vertexSizeMax-vertexSizeMin)

  }else{

    s <- rep(mean(c(vertexSizeMin,vertexSizeMax)),length(s))

  }

  ##------------------------------------------------------------
  ## Node colours
  ##------------------------------------------------------------

  cols <- viridisLite::viridis(100)

  colour <- cols[
    cut(
      V(g)$trophic_level,
      breaks = 100,
      labels = FALSE
    )
  ]

  ##------------------------------------------------------------
  ## Nodes
  ##------------------------------------------------------------

  x_scale <- if(modules) module_spacing else module_spacing * 1.5

  nodes <- data.frame(

    id = seq_len(vcount(g)),

    label = labels,

    title = paste0(
      "<b>",V(g)$name,"</b><br>",
      "Trophic level: ",round(V(g)$trophic_level,2),"<br>",
      "Degree: ",degree(g),"<br>",
      "Module: ",V(g)$module
    ),

    size = s,

    color.background = colour,
    color.border = "black",

    x = x*x_scale,

    y = y*180,

    stringsAsFactors = FALSE

  )

  ##------------------------------------------------------------
  ## Physics mode
  ##------------------------------------------------------------
  ## "x"    -> los nodos pueden reacomodarse horizontalmente pero
  ##           mantienen fija su posición vertical (trophic level)
  ## "full" -> física completa en x e y, con anclaje central y
  ##           auto-apagado post-estabilización para evitar el giro
  ##           continuo característico de forceAtlas2Based

  if(physics=="x"){

    nodes$physics <- TRUE
    nodes$fixed.x <- FALSE
    nodes$fixed.y <- TRUE

  }else{

    nodes$physics <- TRUE
    nodes$fixed.x <- FALSE
    nodes$fixed.y <- FALSE

  }

  ##------------------------------------------------------------
  ## Edges
  ##------------------------------------------------------------

  edge_list <- ends(g,E(g),names=FALSE)

  edges <- data.frame(

    from = edge_list[,1],
    to   = edge_list[,2],

    arrows = "to",

    smooth = FALSE,

    color = "gray70",

    stringsAsFactors = FALSE

  )

  ##------------------------------------------------------------
  ## Plot
  ##------------------------------------------------------------

  p <- visNetwork(
    nodes,
    edges,
    width = "100%",
    height = "850px"
  ) |>

    visNodes(
      shape = "dot",
      borderWidth = 1,
      font = list(
        size = label_size,
        face = "Arial"
      )
    ) |>

    visEdges(
      arrows = list(
        to = list(enabled = TRUE, scaleFactor = 0.4)
      ),
      color = list(color = "gray70"),
      smooth = FALSE
    ) |>

    visOptions(
      highlightNearest = list(
        enabled = highlight,
        degree = 1,
        hover = TRUE
      ),
      nodesIdSelection = search
    ) |>

    visInteraction(
      dragNodes = TRUE,
      dragView = TRUE,
      zoomView = TRUE,
      navigationButtons = TRUE
    ) |>

    visPhysics(

      enabled = TRUE,

      solver = "forceAtlas2Based",

      stabilization = list(

        enabled = TRUE,

        iterations = 1000,

        updateInterval = 25

      ),

      forceAtlas2Based = list(

        gravitationalConstant = -50,

        centralGravity = 0.01,

        springLength = 120,

        springConstant = 0.05,

        damping = 0.6,

        avoidOverlap = 0.5

      )

    )

  p <- p |>
    visEvents(
      stabilizationIterationsDone = "function () {
        this.setOptions({ physics: false });
      }"
    )

  return(p)

}
