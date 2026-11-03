# Generate directed Erdos-Renyi random networks with at least 1 basal node and only one component

This uses the igraph's function sample_gnm to generate nsim random
networks with the same number of nodes and links than the parameter ig
and two restrictions:

## Usage

``` r
generate_er_basal(ig, nsim = 1000)
```

## Arguments

- ig:

  igraph object with parameters to use in the random network
  simulations: number of species/nodes and number of links/edges

- nsim:

  number of simulations

## Value

a list with igraph objects

## Details

1.  at least one basal species/node, that is a species that has no prey.

2.  only 1 connected component so there is no disconnected species or
    sub-community.

## Examples

``` r

generateERbasal(netData[[1]])
```
