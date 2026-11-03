# 3D Visualization of a Multiplex Ecological Network

This function generates a 3D visualization of a multiplex network
composed of multiple layers (e.g., trophic, mutualistic, and competitive
networks) using `rgl`. Layers are displayed stacked along the Z-axis.
Node positions are shared across layers, and nodes can be colored by
module membership or functional type. An optional aggregate layer
(merging all layers) can also be displayed.

## Usage

``` r
plot_multi3D(
  g.list,
  layer.colors,
  as.undirected = T,
  layer.layout = "auto",
  layer.labels = "auto",
  layer.labels.cex = 2,
  edge.colors = "auto",
  edge.normalize = F,
  edge.size.scale = 1,
  node.colors = "auto",
  node.size.values = 0.5,
  node.size.scale = 1,
  node.alpha = 1,
  edge.alpha = 1,
  layer.alpha = "auto",
  layout = "fr",
  show.nodeLabels = F,
  show.aggregate = F,
  aggr.alpha = "auto",
  aggr.color = "#dadada",
  node.colors.aggr = "#dadada",
  layer.scale = 2,
  layer.shift.x = 0,
  layer.shift.y = 0,
  layer.space = 1.5,
  FOV = 30
)
```

## Arguments

- g.list:

  A list of `igraph` objects, each representing a layer of the multiplex
  network.

- layer.colors:

  A vector of colors corresponding to each network layer.

- as.undirected:

  Logical; if `TRUE`, converts each layer to undirected for
  visualization purposes (default: `TRUE`).

- layer.layout:

  A layout matrix shared by all layers, or `"auto"` to use
  `layoutMultiplex()` (default: `"auto"`).

- layer.labels:

  A character vector for labeling layers. `"auto"` assigns "Layer 1",
  "Layer 2", etc. (default: `"auto"`).

- layer.labels.cex:

  Numeric; text size for layer labels (default: `2`).

- edge.colors:

  Edge colors (a single color or vector per layer), or `"auto"` to match
  layer color (default: `"auto"`).

- edge.normalize:

  Logical; if `TRUE`, normalizes edge weights using logarithmic scale
  (default: `FALSE`).

- edge.size.scale:

  Numeric or vector to control edge thickness (default: `1`).

- node.colors:

  `"auto"` to use layer colors, or an `N x L` matrix for node colors per
  layer.

- node.size.values:

  `"auto"` to size nodes by strength, or a fixed numeric value (default:
  `0.5`).

- node.size.scale:

  Numeric or vector to scale node size (default: `1`).

- node.alpha:

  Transparency for nodes (0–1; default: `1`).

- edge.alpha:

  Transparency for edges (0–1; default: `1`).

- layer.alpha:

  A numeric vector with layer transparency, or `"auto"` to set 0.5 for
  all layers.

- layout:

  A character string indicating the layout method for
  `layoutMultiplex()` (default: `"fr"`).

- show.nodeLabels:

  Logical; if `TRUE`, node labels are shown (default: `FALSE`).

- show.aggregate:

  Logical; if `TRUE`, a layer with the aggregated network is shown
  (default: `FALSE`).

- aggr.alpha:

  Transparency for the aggregate layer (default: `"auto"` which is set
  to 0.5).

- aggr.color:

  Color of the background layer for the aggregate network (default:
  `"#dadada"`).

- node.colors.aggr:

  Node color used in the aggregate network (default: `"#dadada"`).

- layer.scale:

  Scale factor for the extent of each layer (default: `2`).

- layer.shift.x:

  Horizontal X shift per layer (default: `0`).

- layer.shift.y:

  Horizontal Y shift per layer (default: `0`).

- layer.space:

  Vertical spacing between layers along Z (default: `1.5`).

- FOV:

  Field of view for the `rgl` perspective visualization (default: `30`).

## Details

Each layer is visualized as a semi-transparent plane, on which the
respective network is plotted. Nodes are shared across layers and
maintain consistent positions. Optionally, the aggregated network (union
of all layers) can be visualized at the top. This function is designed
for ecological multiplex networks and facilitates interpretation of
modular structure, interaction types, or perturbation effects.

This function is adapted from the `plot_multiplex3D()` function in the
[`muxViz`](https://manlius.github.io/muxViz/) project.

## Examples

``` r
if (FALSE) { # \dontrun{
library(multiweb)
fileName <- system.file("extdata", package = "multiweb")
dn <- list.files(fileName, pattern = "^Kefi2015.*\\.txt$")
g_list <- readNetwork(dn, fileName, skipColumn = 2)
names(g_list) <- c("Negative", "Positive", "Trophic")


# Visualize with community-colored nodes
plot_multi3D(
  g.list = g_list,
  layer.colors = RColorBrewer::brewer.pal(3, "Set2"),
  node.colors = "auto",
  node.size.scale = 0.5,
  node.size.values = "auto",
  node.colors.aggr = "#999999",
  show.aggregate = TRUE
)
} # }
```
