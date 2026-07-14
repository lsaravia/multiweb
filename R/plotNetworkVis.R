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
plot_troph_level_visNet <- function(
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


#' Plot a Multiplex Ecological Network Organized by Trophic Level using visNetwork
#'
#' Extends [plot_troph_level_visNet()] to multiplex ecological networks, where
#' species interact through several distinct types of interactions (e.g. trophic vs.
#' non-trophic: facilitation, competition, commensalism, amensalism) represented as
#' separate "layers" sharing the same set of nodes. Trophic level is computed using
#' only the trophic layer, while all layers are drawn together with distinguishable
#' edge colors/styles and an interactive layer filter.
#'
#' @param g Either (a) a single `igraph` object with an edge attribute (named by
#'   `layer_attr`) identifying the interaction type of each edge, or (b) a named list
#'   of `igraph` objects, one per layer, all sharing (a subset of) the same vertex set
#'   via `V(g)$name`. List names are used as layer labels.
#' @param layer_attr Name of the edge attribute holding the layer/interaction-type
#'   label, used only when `g` is a single `igraph` object (default: `"layer"`).
#' @param trophic_layer Value identifying which layer represents trophic (feeding)
#'   interactions; used to compute trophic level via `NetIndices::TrophInd()`
#'   (default: `"trophic"`). Species with no trophic-layer edges (e.g. species that
#'   only participate in non-trophic interactions) are assigned trophic level = 1
#'   for display purposes, with a warning.
#' @param vertexSizeMin,vertexSizeMax Minimum/maximum node size (defaults: 8, 25).
#' @param modules Logical; if `TRUE`, nodes are grouped by community modules along
#'   the x-axis, detected on the flattened (all-layers-pooled) network (default: `FALSE`).
#' @param weights Edge weights used for community detection when `modules = TRUE`
#'   and `community_obj` is not supplied (default: `NA`).
#' @param community_obj Optional pre-computed community detection object.
#' @param node_weights Optional numeric vector used to size nodes. If `NA` (default),
#'   total node degree pooled across all layers is used, so multiplex hubs are sized
#'   according to their combined role across interaction types.
#' @param module_spacing Horizontal spacing multiplier between community modules
#'   (default: 350).
#' @param physics Layout mode, `"x"` (default) or `"full"` — see
#'   [plot_troph_level_visNet()] for details.
#' @param highlight Logical; hover-highlight nearest neighbors (default: `TRUE`).
#' @param search Logical; add node-search dropdown (default: `TRUE`).
#' @param label_size Node label font size in pixels (default: 16).
#' @param layer_colors Optional named character vector mapping layer names to edge
#'   colors. If `NULL` (default), a `viridisLite` palette is generated automatically.
#' @param layer_dashed Optional character vector of layer names whose edges should be
#'   drawn dashed, to distinguish them from solid edges regardless of color (typically
#'   used to mark non-trophic layers, e.g. `layer_dashed = c("facilitation","competition")`).
#' @param layer_filter Logical; if `TRUE` (default), adds an interactive dropdown to
#'   filter/isolate edges by layer.
#'
#' @return A `visNetwork` HTML widget.
#'
#' @export
#'
#' @importFrom igraph as_adjacency_matrix V E "E<-" degree cluster_spinglass membership
#'   vertex_attr_names edge_attr_names edge_attr vcount ecount is_directed is_igraph
#'   graph_from_data_frame as_data_frame subgraph_from_edges
#' @importFrom NetIndices TrophInd
#' @importFrom visNetwork visNetwork visNodes visEdges visOptions visInteraction
#'   visPhysics visEvents visLegend
#' @importFrom viridisLite viridis
plot_troph_level_visNet_multi <- function(
    g,
    layer_attr = "layer",
    trophic_layer = "trophic",
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
    label_size = 16,
    layer_colors = NULL,
    layer_dashed = NULL,
    layer_filter = TRUE
){

  physics <- match.arg(physics)

  ##------------------------------------------------------------
  ## Normalize input: single igraph + layer edge attribute
  ##------------------------------------------------------------
  ## Accepts either a single igraph with a `layer_attr` edge attribute, or a named
  ## list of igraph objects (one per layer) sharing the same vertex names. Both are
  ## flattened into one igraph (`g_full`) with a `layer` edge attribute, which is
  ## what the rest of the function operates on.

  if(is.list(g) && !is_igraph(g)){

    layer_names <- names(g)

    if(is.null(layer_names) || any(layer_names==""))
      stop("If `g` is a list of graphs, it must be a fully named list (names = layer labels).")

    all_node_names <- unique(unlist(lapply(g, function(x) V(x)$name)))

    directed <- any(vapply(g, is_directed, logical(1)))

    el_list <- lapply(layer_names, function(ln){

      el <- as_data_frame(g[[ln]], what = "edges")

      if(nrow(el) > 0) el$layer <- ln

      el[, c("from","to","layer")]

    })

    edges_all <- do.call(rbind, el_list)

    g_full <- graph_from_data_frame(
      edges_all,
      directed = directed,
      vertices = data.frame(name = all_node_names, stringsAsFactors = FALSE)
    )

  }else{

    if(!layer_attr %in% edge_attr_names(g))
      stop(sprintf(
        "Edge attribute '%s' not found in `g`. Set `layer_attr` to the correct attribute name.",
        layer_attr
      ))

    g_full <- g
    E(g_full)$layer <- edge_attr(g_full, layer_attr)

  }

  layer_names_all <- sort(unique(E(g_full)$layer))

  trophic_match <- which(tolower(layer_names_all) == tolower(trophic_layer))

  if(length(trophic_match) == 0)
    stop(sprintf(
      "`trophic_layer` = '%s' not found among layers: %s",
      trophic_layer, paste(layer_names_all, collapse=", ")
    ))

  trophic_layer <- layer_names_all[trophic_match[1]]

  ##------------------------------------------------------------
  ## Trophic level (computed from the trophic layer only)
  ##------------------------------------------------------------

  trophic_eids <- which(E(g_full)$layer == trophic_layer)

  g_troph <- subgraph_from_edges(g_full, eids = trophic_eids, delete.vertices = FALSE)

  adj <- as_adjacency_matrix(g_troph, sparse = FALSE)
  tl <- NetIndices::TrophInd(adj)

  troph_level <- tl$TL

  if(any(!is.finite(troph_level))){

    warning(
      "Some species have no trophic-layer edges (e.g. purely non-trophic taxa); ",
      "assigning trophic level = 1 for display purposes."
    )

    troph_level[!is.finite(troph_level)] <- 1

  }

  V(g_full)$trophic_level <- troph_level

  ##------------------------------------------------------------
  ## Node weights (pooled multiplex degree by default)
  ##------------------------------------------------------------

  if(length(node_weights)==1 && is.na(node_weights)){

    node_weights <- degree(g_full, mode="all")

  }

  V(g_full)$weight <- node_weights

  ##------------------------------------------------------------
  ## Community modules (detected on the flattened multiplex network)
  ##------------------------------------------------------------

  if(modules){

    if(is.null(community_obj))
      community_obj <- cluster_spinglass(g_full, weights = weights)

    V(g_full)$module <- membership(community_obj)

  }else{

    V(g_full)$module <- 1

  }

  ##------------------------------------------------------------
  ## Labels
  ##------------------------------------------------------------

  labels <-
    if("label" %in% vertex_attr_names(g_full))
      V(g_full)$label
  else
    V(g_full)$name

  ##------------------------------------------------------------
  ## Coordinates
  ##------------------------------------------------------------

  x <- if(modules){

    jitter(as.numeric(V(g_full)$module), amount = .45)

  }else{

    runif(vcount(g_full), min = -1, max = 1)

  }

  y <- -V(g_full)$trophic_level

  ##------------------------------------------------------------
  ## Node size
  ##------------------------------------------------------------

  s <- sqrt(V(g_full)$weight)

  if(diff(range(s)) > 0){

    s <- (s-min(s))/(max(s)-min(s))

    s <- vertexSizeMin + s*(vertexSizeMax-vertexSizeMin)

  }else{

    s <- rep(mean(c(vertexSizeMin,vertexSizeMax)),length(s))

  }

  ##------------------------------------------------------------
  ## Node colours (by trophic level, as in the base function)
  ##------------------------------------------------------------

  cols <- viridisLite::viridis(100)

  colour <- cols[
    cut(V(g_full)$trophic_level, breaks = 100, labels = FALSE)
  ]

  ##------------------------------------------------------------
  ## Nodes
  ##------------------------------------------------------------

  x_scale <- if(modules) module_spacing else module_spacing * 1.5

  nodes <- data.frame(

    id = seq_len(vcount(g_full)),

    label = labels,

    title = paste0(
      "<b>",V(g_full)$name,"</b><br>",
      "Trophic level: ",round(V(g_full)$trophic_level,2),"<br>",
      "Multiplex degree: ",degree(g_full),"<br>",
      "Module: ",V(g_full)$module
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
  ## Layer colors / line styles
  ##------------------------------------------------------------

  if(is.null(layer_colors)){

    pal <- viridisLite::viridis(length(layer_names_all), option = "turbo")
    layer_colors <- setNames(pal, layer_names_all)

  }else{

    missing_layers <- setdiff(layer_names_all, names(layer_colors))

    if(length(missing_layers) > 0)
      stop("`layer_colors` is missing entries for: ", paste(missing_layers, collapse=", "))

  }

  if(is.null(layer_dashed)) layer_dashed <- character(0)

  ##------------------------------------------------------------
  ## Edges
  ##------------------------------------------------------------
  ## Each layer gets its own color, an optional dashed style, and a slightly
  ## different curvature (based on layer order) so that overlapping edges between
  ## the same pair of species across different layers remain visually separable.

  edge_df <- as_data_frame(g_full, what = "edges")

  id_map <- setNames(seq_len(vcount(g_full)), V(g_full)$name)

  layer_idx <- match(edge_df$layer, layer_names_all)

  edges <- data.frame(

    from = id_map[edge_df$from],
    to   = id_map[edge_df$to],

    arrows = "to",

    color = layer_colors[edge_df$layer],
    dashes = edge_df$layer %in% layer_dashed,

    title = paste0("<b>",edge_df$layer,"</b>"),
    group = edge_df$layer,

    smooth.enabled = TRUE,
    smooth.type = "curvedCW",
    smooth.roundness = 0.1 + 0.15*(layer_idx-1),

    stringsAsFactors = FALSE

  )

  ##------------------------------------------------------------
  ## Legend (one entry per layer)
  ##------------------------------------------------------------

  legend_edges <- data.frame(
    label = layer_names_all,
    color = unname(layer_colors[layer_names_all]),
    dashes = layer_names_all %in% layer_dashed,
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
      smooth = FALSE  # per-edge smooth.* columns above take precedence
    ) |>

    visOptions(
      highlightNearest = list(
        enabled = highlight,
        degree = 1,
        hover = TRUE
      ),
      nodesIdSelection = search,
      selectedBy = if(layer_filter) list(variable = "group", multiple = TRUE) else NULL
    ) |>

    visLegend(
      addEdges = legend_edges,
      useGroups = FALSE,
      position = "right"
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

  if(physics == "full"){

    p <- p |>
      visEvents(
        stabilizationIterationsDone = "function () {
          this.setOptions({ physics: false });
        }"
      )

  }

  return(p)

}
