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
