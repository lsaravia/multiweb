# Singular Value Decomposition (SVD) Entropy and Rank Calculation

This function calculates the entropy and rank of a network using
Singular Value Decomposition (SVD). It accepts either an adjacency
matrix or an igraph object and processes it accordingly.

## Usage

``` r
calc_svd_entropy(A, threshold_factor = 1e-06)
```

## Arguments

- A:

  An adjacency matrix (square, weighted or unweighted) or an igraph
  object.

- threshold_factor:

  A small value (default = 1e-6) to determine numerical rank
  approximation.

## Value

A list containing:

- Rank:

  Effective rank of the matrix.

- Entropy:

  SVD entropy of the network.
