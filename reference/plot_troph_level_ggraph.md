# Plot an Ecological Network Organized by Trophic Level using ggraph

This function visualizes a food web using `ggplot2`, where nodes
represent species and edges represent interactions. Nodes are positioned
by trophic level on the y-axis, and optionally, community modules on the
x-axis. Is similar to
[`plot_troph_level_ggplot()`](https://lsaravia.github.io/multiweb/reference/plot_troph_level_ggplot.md),
but uses `ggraph` for layout and rendering, which make better arrows.

## Usage

``` r
plot_troph_level_ggraph(
  g,
  vertexSizeFactor = 5,
  vertexSizeMin = 5,
  modules = FALSE,
  weights = NA,
  community_obj = NULL,
  use_numbers = FALSE,
  label_size = 4,
  arrow_size = 0.1,
  node_weights = NA
)
```

## Arguments

- g:

  An `igraph` network.

- vertexSizeFactor:

  Factor for node size scaling (default 5).

- vertexSizeMin:

  Minimum node size (default 5).

- modules:

  If TRUE, uses `cluster_spinglass` for community detection.

- weights:

  Edge weights for clustering (default NA).

- community_obj:

  Optional community detection object.

- use_numbers:

  Label nodes by numeric ID.

- label_size:

  Text size.

- arrow_size:

  Arrowhead size.

- node_weights:

  Optional numeric vector for node weights. If NA, uses degree. If it is
  NULL and the input graph has a ‘weight’ vertex attribute, then that
  attribute will be used.

## Value

A `ggraph` plot.
