# Generate a Sequence of Shuffled Networks and Track SVD/Modularity Evolution

This function generates a sequence of shuffled networks from an original
graph by applying incremental shuffling steps. It tracks modularity and
SVD entropy throughout the process, and stops when the change in these
metrics falls below a specified `tolerance` threshold.

## Usage

``` r
generate_shuffled_seq_tol(
  input_graph,
  delta = 10,
  max_iterations = 100,
  tolerance = 0.001,
  directed = TRUE,
  weighted = TRUE,
  modularity_func = cluster_infomap,
  shuffle_func = shuffle_network_deg
)
```

## Arguments

- input_graph:

  An igraph object or adjacency matrix (directed or undirected, weighted
  or unweighted).

- delta:

  Integer. Number of link swaps per iteration (default: 10).

- max_iterations:

  Integer. Maximum number of iterations until network metrics stabilize
  (default: 100).

- tolerance:

  Numeric. Tolerance for metric stability - minimum change in modularity
  or SVD entropy (default: 1e-3).

- directed:

  Logical. Whether the network is directed (default: TRUE).

- weighted:

  Logical. Whether the network is weighted (default: TRUE).

- modularity_func:

  Function. Function to calculate modularity (default:
  `cluster_infomap`).

- shuffle_func:

  Function. Function to perform network shuffling (default:
  `shuffle_network_deg`).

## Value

A list containing:

- New_A:

  Shuffled adjacency matrix

- Metrics:

  Data frame tracking metric evolution across iterations, including the
  original network
