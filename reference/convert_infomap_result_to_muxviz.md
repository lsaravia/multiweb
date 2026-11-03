# Convert Infomap Community Output to muxViz-Compatible Format

This function converts the community detection result from
[`run_infomap_multi()`](https://lsaravia.github.io/multiweb/reference/run_infomap_multi.md)
(from the `multiweb` package) into a format compatible with `muxViz`,
allowing direct use with visualization function `plot_multimodules()`.

## Usage

``` r
convert_infomap_result_to_muxviz(communities)
```

## Arguments

- communities:

  A data frame returned by `run_infomap_multi()$communities`, containing
  columns `node`, `layer`, `module`, and `flow`.

## Value

A named list with the following components:

- membership.multi:

  A data frame with columns `node`, `layer`, and `module` representing
  the state node assignment. Node and layer are returned as integers
  (1-based index).

- membership.aggr:

  A data frame with columns `node` and `module`, where each physical
  node is assigned to its dominant module across all layers.

- modules.multi:

  Integer; the number of unique modules detected in the multilayer
  network.

- modules.aggr:

  Integer; the number of unique modules in the aggregated network.

## Details

It creates two key data frames:

- `membership.multi`: Community membership for each node-layer
  combination.

- `membership.aggr`: Aggregate community membership per physical node
  (across layers), based on the highest flow.

## See also

[`run_infomap_multi`](https://lsaravia.github.io/multiweb/reference/run_infomap_multi.md),
[`plot_multi3D`](https://lsaravia.github.io/multiweb/reference/plot_multi3D.md),
`plot_multimodules`

## Examples

``` r
if (FALSE) { # \dontrun{
# Run Infomap community detection on a multiplex ecological network
g_list <- readNetwork(dn, fileName, skipColumn = 2)
names(g_list) <- c("Negative", "Positive", "Trophic")
result <- run_infomap_multi(g_list, layer_names = names(g_list))

# Convert for muxViz plotting
mux_data <- convert_infomap_result_to_muxviz(result$communities)

# Visualize
plot_multimodules(mux_data)
} # }
```
