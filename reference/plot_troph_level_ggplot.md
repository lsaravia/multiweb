# Plot an Ecological Network Organized by Trophic Level

This function visualizes a food web using `ggplot2`, where nodes
represent species and edges represent interactions. Nodes are positioned
by trophic level on the y-axis, and optionally, community modules on the
x-axis.

## Usage

``` r
plot_troph_level_ggplot(
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
)
```

## Arguments

- g:

  An `igraph` object representing the ecological network.

- vertexSizeFactor:

  Numeric factor to scale node size (default: 5).

- vertexSizeMin:

  Minimum node size (default: 5).

- modules:

  Logical; if `TRUE`, nodes are grouped by community modules (default:
  `FALSE`).

- weights:

  Edge weights for community detection (default: `NA`).

- node_weights:

  Optional numeric vector for node weights. If `NULL` and the graph has
  a vertex attribute called `weight`, that attribute is used. If not,
  node degree is used as fallback.

- community_obj:

  Optional community detection object.

- use_numbers:

  Logical; if `TRUE`, label nodes by numeric ID (default: `FALSE`).

- label_size:

  Label font size (default: 4).

- arrow_size:

  Arrowhead size (default: 0.15).

- shorten_factor:

  Factor to shorten edges for arrows (default: 0.005).

## Value

A `ggplot` object visualizing the trophic structure of the network.
