# Aggregate a Multiplex Network into a Single Layer

Combines a list of `igraph` layers into one aggregate network. Edges
present in multiple layers will have their weights summed.

## Usage

``` r
aggregate_multiplex_network(g_list, directed = TRUE)
```

## Arguments

- g_list:

  A list of `igraph` objects (directed or undirected).

- directed:

  Logical; whether to treat the result as directed (default: `TRUE`).

## Value

An `igraph` object representing the aggregate network.
