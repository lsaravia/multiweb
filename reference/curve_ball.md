# Curve ball algorithm to generate random networks with given an igraph network

This extract the adjacency matrix from an igraph object and generates a
randomized version of the network with the same row and column totals,
the results have the same in and out degree sequence than the original
network. The diagonals are not avoided so it can generate self-links or
cannibalism in the context of food-webs. In the case that the network
has multiple components (or disconnected networks) the algorithm
simulates around the same number of components, if the original network
has 1 component the algorithm enforces that the results have all 1
component. If the edge attribute weight is present and istrength=TRUE
then the weigth is additionally randomized keeping the column sum equal
to the original matrix.

## Usage

``` r
curve_ball(g, nsim = 1000, istrength = FALSE, fix_loops = FALSE)
```

## Arguments

- g:

  igraph object to extract adjacency matrix

- nsim:

  number of generated random networks

- istrength:

  if TRUE the edge attribute weight is taken as interaction strength and
  randomized keeping the column sum equal to the original matrix.

- fix_loops:

  if TRUE the self-loops (the diagonal) are kept fixed: they are not
  traded and no new loops are created. Default FALSE, the diagonal is
  randomized like any other cell. See
  [`shuffle_network_deg()`](https://lsaravia.github.io/multiweb/reference/shuffle_network_deg.md)
  for the motivation.

## Value

a list of randomized igraph objects

## Details

A warning is given if the adjacency matrix has values greater than 1 (a
multigraph), because only the cells equal to 1 are used.

For each randomized network `20 * n` trades (n = number of rows) are
attempted; in the networks tested this is enough for the fraction of
links shared with the original network to stop decreasing, giving the
same distribution as
[`shuffle_network_deg()`](https://lsaravia.github.io/multiweb/reference/shuffle_network_deg.md)
after it converges.

Based on:

## References

Strona, G. et al. 2014. A fast and unbiased procedure to randomize
ecological binary matrices with fixed row and column totals. -Nat. Comm.
5: 4114. doi: 10.1038/ncomms5114

## Examples

``` r
g <- netData[[1]]
nets <- curve_ball(g, nsim = 5)
length(nets)
#> [1] 5

# The randomized networks have the same in and out degree than the original
all(igraph::degree(nets[[1]], mode = "out") == igraph::degree(g, mode = "out"))
#> [1] TRUE
all(igraph::degree(nets[[1]], mode = "in") == igraph::degree(g, mode = "in"))
#> [1] TRUE
```
