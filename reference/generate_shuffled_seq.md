# Generate a Sequence of Shuffled Networks and Track Metric Evolution

This function generates a sequence of shuffled networks from an original
graph, applying incremental shuffling steps while tracking modularity,
SVD rank and entropy .

## Usage

``` r
generate_shuffled_seq(
  original_graph,
  max_delta = 10,
  delta = 10,
  directed = TRUE,
  weighted = TRUE,
  shuffle_func = shuffle_network_ws,
  modularity_func = cluster_infomap
)
```

## Arguments

- original_graph:

  An igraph object representing the original network.

- max_delta:

  Integer. Number of total shuffling steps to perform.

- delta:

  Integer. Number of swaps per shuffling step.

- directed:

  Logical. Whether the network is directed (default: TRUE).

- weighted:

  Logical. Whether the network is weighted (default: TRUE).

- shuffle_func:

  Function. A network shuffling function (default:
  `shuffle_network_ws`).

- modularity_func:

  Function. A modularity calculation function (default:
  `cluster_infomap`).

## Value

A list with:

- Networks:

  A list of igraph objects representing the shuffled networks at each
  step.

- Metrics:

  A tibble containing the step number, SVD entropy, SVD Rank, and
  modularity score.

## Details

`shuffle_func` is called as
`shuffle_func(A, delta = delta, directed = directed, weighted = weighted)`
where `A` is the adjacency matrix of the current network, so a wrapper
must accept those arguments.

## Examples

``` r
g <- netData[[1]]
# Keep the self-loops (cannibalism) fixed during the randomization
res <- generate_shuffled_seq(g, max_delta = 60, delta = 10,
  shuffle_func = function(x, delta, directed, weighted)
    shuffle_network_deg(x, delta, directed, weighted, fix_loops = TRUE))
head(res$Metrics)
#> # A tibble: 6 × 4
#>    Step SVD_Entropy SVD_Rank Modularity
#>   <dbl>       <dbl>    <int>      <dbl>
#> 1     1       0.886       21          0
#> 2     2       0.890       21          0
#> 3     3       0.900       22          0
#> 4     4       0.903       22          0
#> 5     5       0.901       22          0
#> 6     6       0.900       22          0
```
