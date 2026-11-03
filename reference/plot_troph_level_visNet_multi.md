# Plot a Multiplex Ecological Network Organized by Trophic Level using visNetwork

Extends
[`plot_troph_level_visNet()`](https://lsaravia.github.io/multiweb/reference/plot_troph_level_visNet.md)
to multiplex ecological networks, where species interact through several
distinct types of interactions (e.g. trophic vs. non-trophic:
facilitation, competition, commensalism, amensalism) represented as
separate "layers" sharing the same set of nodes. Trophic level is
computed using only the trophic layer, while all layers are drawn
together with distinguishable edge colors/styles and an interactive
layer filter.

## Usage

``` r
plot_troph_level_visNet_multi(
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
  physics = c("x", "full"),
  highlight = TRUE,
  search = TRUE,
  label_size = 16,
  layer_colors = NULL,
  layer_dashed = NULL,
  layer_filter = TRUE,
  legend = TRUE,
  width = "100%",
  height = "850px"
)
```

## Arguments

- g:

  Either (a) a single `igraph` object with an edge attribute (named by
  `layer_attr`) identifying the interaction type of each edge, or (b) a
  named list of `igraph` objects, one per layer, all sharing (a subset
  of) the same vertex set via `V(g)$name`. List names are used as layer
  labels.

- layer_attr:

  Name of the edge attribute holding the layer/interaction-type label,
  used only when `g` is a single `igraph` object (default: `"layer"`).

- trophic_layer:

  Value identifying which layer represents trophic (feeding)
  interactions; used to compute trophic level via
  [`NetIndices::TrophInd()`](https://rdrr.io/pkg/NetIndices/man/TrophInd.html)
  (default: `"trophic"`). Species with no trophic-layer edges (e.g.
  species that only participate in non-trophic interactions) are
  assigned trophic level = 1 for display purposes, with a warning.

- vertexSizeMin, vertexSizeMax:

  Minimum/maximum node size (defaults: 8, 25).

- modules:

  Logical; if `TRUE`, nodes are grouped by community modules along the
  x-axis, detected on the flattened (all-layers-pooled) network
  (default: `FALSE`).

- weights:

  Edge weights used for community detection when `modules = TRUE` and
  `community_obj` is not supplied (default: `NA`).

- community_obj:

  Optional pre-computed community detection object.

- node_weights:

  Optional numeric vector used to size nodes. If `NA` (default), total
  node degree pooled across all layers is used, so multiplex hubs are
  sized according to their combined role across interaction types.

- module_spacing:

  Horizontal spacing multiplier between community modules (default:
  350).

- physics:

  Layout mode, `"x"` (default) or `"full"` — see
  [`plot_troph_level_visNet()`](https://lsaravia.github.io/multiweb/reference/plot_troph_level_visNet.md)
  for details. In both modes, physics is disabled and node positions are
  unlocked once stabilization finishes, so the trophic-level ordering
  (`"x"`) or force-directed layout (`"full"`) is only the starting
  point; afterwards any node can be freely dragged in any direction.

- highlight:

  Logical; hover-highlight nearest neighbors (default: `TRUE`).

- search:

  Logical; add node-search dropdown (default: `TRUE`).

- label_size:

  Node label font size in pixels (default: 16).

- layer_colors:

  Optional named character vector mapping layer names to edge colors. If
  `NULL` (default), a `viridisLite` palette is generated automatically.

- layer_dashed:

  Optional character vector of layer names whose edges should be drawn
  dashed, to distinguish them from solid edges regardless of color
  (typically used to mark non-trophic layers, e.g.
  `layer_dashed = c("facilitation","competition")`).

- layer_filter:

  Logical; if `TRUE` (default), adds an interactive dropdown to
  filter/isolate edges by layer.

- legend:

  Logical; if `TRUE` (default), displays the layer-type legend.

- width:

  Width of the `visNetwork` widget (default: `"100%"`). Accepts any CSS
  width string.

- height:

  Height of the `visNetwork` widget (default: `"850px"`). Accepts any
  CSS height string. When embedding in slide frameworks like `xaringan`,
  set an explicit pixel height small enough to fit the slide (e.g.
  `"400px"`) — the widget does not shrink to fit its container
  automatically.

## Value

A `visNetwork` HTML widget.

## Details

When present, existing `igraph` visual attributes are respected: vertex
`color` and `size`, and edge `color` and `width`.
