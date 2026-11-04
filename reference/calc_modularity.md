# Calculation of Modularity for a List of igraph Objects

This function calculates the modularity of a list of networks provided
in the `ig` parameter. Modularity measures the strength of division of a
network into modules (communities).

## Usage

``` r
calc_modularity(
  ig,
  ncores = 0,
  cluster_function = cluster_spinglass,
  directed = FALSE
)
```

## Arguments

- ig:

  A list of igraph objects for which modularity will be calculated.

- ncores:

  Number of cores used for parallel computation. If set to 0, sequential
  processing is used.

- cluster_function:

  The **igraph** clustering function to use, or a function that returns
  an igraph communities object (e.g.
  [`run_infomap()`](https://lsaravia.github.io/multiweb/reference/run_infomap.md)).
  Default `cluster_spinglass`.

- directed:

  Logical; if `TRUE` and the network is directed, the directed
  modularity of the partition is computed. Default `FALSE` (undirected
  modularity), which is consistent with the default
  `cluster_spinglass()` because it ignores the direction of the edges.
  Use `TRUE` only with a clustering function that takes the direction
  into account, e.g. `run_infomap(directed = TRUE)`.

## Value

A data frame with a column `Modularity` containing the modularity values
for each network.

## Details

**Note:** `cluster_spinglass()` only works on networks with a single
connected component.

The modularity is always recomputed from the membership of the partition
with
[`igraph::modularity()`](https://r.igraph.org/reference/modularity.igraph.html),
and not taken from the value stored in the object returned by the
clustering function, because that value is not comparable between
algorithms (for `cluster_spinglass()` it is not the weighted modularity
of the partition). If a network has a 'weight' edge attribute it is used
both by the clustering function and to compute the (weighted)
modularity.

## Examples

``` r
if (FALSE) { # \dontrun{
nullg <- generateERbasal(netData[[1]], 10)
calc_modularity(nullg)

# Directed modularity with an algorithm that uses the direction of the edges
calc_modularity(netData[[1]], cluster_function = run_infomap, directed = TRUE)
} # }
```
