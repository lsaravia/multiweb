# Random Network Rewiring Without Preserving Degree Distribution

This function randomly rewires a directed network while preserving the
total number of links but NOT the degree distribution, this procedure is
adapted from Huaylla et al. (2024) and is based on an approach described
by Watts and Strogatz (1998) for small-world networks.

## Usage

``` r
shuffle_network_ws(input_graph, delta = 100, directed = TRUE, weighted = TRUE)
```

## Arguments

- input_graph:

  A square adjacency matrix (directed, weighted or unweighted) or an
  igraph object.

- delta:

  Number of rewiring attempts to perform.

- directed:

  Logical, whether the network is directed (default: TRUE).

- weighted:

  Logical, whether the network is weighted (default: TRUE).

## Value

A shuffled adjacency matrix preserving the number of links but not
degree distribution.

## Details

The procedure is: we select a single pair of connected nodes and
disconnect one of the nodes at one end of the link. We then connect the
free end to another node that is randomly selected.

## References

Watts, D. J., & Strogatz, S. H. (1998). Collective dynamics of
'small-world' networks. Nature, 393(6684), 440-442.
[doi:10.1038/30918](https://doi.org/10.1038/30918)

Huaylla, C. A., Kuperman, M. N., & Garibaldi, L. A. (2024). Comparison
of two statistical measures of complexity applied to ecological
bipartite networks. Physica A: Statistical Mechanics and Its
Applications, 642, 129764. https://doi.org/10.1016/j.physa.2024.129764
