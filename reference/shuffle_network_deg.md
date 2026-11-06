# Degree-Preserving Network Shuffling via Edge Swaps

This function shuffles a directed network while preserving in-degree and
out-degree. It follows a controlled randomization process by swapping
edges iteratively. Similar to the
[`curve_ball()`](https://lsaravia.github.io/multiweb/reference/curve_ball.md)
function but more computationally demanding. This function is useful for
generating increasingly randomized networks, to generate fully
randomized networks, the
[`curve_ball()`](https://lsaravia.github.io/multiweb/reference/curve_ball.md)
function is preferred.

## Usage

``` r
shuffle_network_deg(
  input_graph,
  delta = 100,
  directed = TRUE,
  weighted = TRUE,
  fix_loops = FALSE
)
```

## Arguments

- input_graph:

  An igraph object or a square adjacency matrix (directed, weighted or
  unweighted).

- delta:

  Number of edge swaps to perform.

- directed:

  Logical. Whether the network is directed (default: TRUE).

- weighted:

  Logical. Whether the network is weighted (default: TRUE). When `FALSE`
  a warning is given if the matrix has values greater than 1 (a
  multigraph).

- fix_loops:

  Logical. If `TRUE` the self-loops (the diagonal) are kept fixed during
  the randomization (default: `FALSE`, the diagonal is randomized like
  any other cell).

## Value

A shuffled adjacency matrix preserving in-degree and out-degree.

## Details

The function preserves the total in-degree and out-degree (row and
column totals) of the adjacency matrix, but it does not preserve the
degree of each layer of a multilayer network nor the reciprocity:
reciprocal pairs A-\>B / B-\>A (e.g. mutualism, competition) are
separated by the swaps.

By default the diagonal is not avoided, so the swaps can create and
destroy self-loops (cannibalism in food-webs). In food-webs with 14-15
loops only 5-8 were left after about two swaps per link. A loop always
falls inside its own module, so losing loops lowers the modularity of
the randomized networks and makes the empirical network look more
modular than it is. With `fix_loops = TRUE` the loops are kept fixed:
they are never selected for a swap and a swap that would create a new
loop is rejected.

## References

Huaylla, C. A., Nacif, M. E., Coulin, C., Kuperman, M. N., & Garibaldi,
L. A. (2021). Decoding information in multilayer ecological networks:
The keystone species case. Ecological Modelling, 460, 109734.
https://doi.org/10.1016/j.ecolmodel.2021.109734

Strona, G., Nappo, D., Boccacci, F., Fattorini, S., & San-Miguel-Ayanz,
J. (2014). A fast and unbiased procedure to randomize ecological binary
matrices with fixed row and column totals. Nature Communications, 5,
4114. https://doi.org/10.1038/ncomms5114
