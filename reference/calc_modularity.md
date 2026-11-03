# Calculation of Modularity for a List of igraph Objects

This function calculates the modularity of a list of networks provided
in the `ig` parameter. Modularity measures the strength of division of a
network into modules (communities).

## Usage

``` r
calc_modularity(ig, ncores = 0, cluster_function = cluster_spinglass)
```

## Arguments

- ig:

  A list of igraph objects for which modularity will be calculated.

- ncores:

  Number of cores used for parallel computation. If set to 0, sequential
  processing is used.

- cluster_function:

  The **igraph** clustering function to use, or a function that returns
  modularity. default `cluster_spinglass`.

  If a network has a 'weight' edge attribute it is used both by the
  clustering function and to compute the (weighted) modularity of the
  resulting partition.

## Value

A data frame with a column `Modularity` containing the modularity values
for each network.

## Details

**Note:** `cluster_spinglass()` only works on networks with a single
connected component.

## Examples

``` r
if (FALSE) { # \dontrun{
nullg <- generateERbasal(netData[[1]], 10)
calc_modularity(nullg)
} # }
```
