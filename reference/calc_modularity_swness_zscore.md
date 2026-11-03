# Calculation of Modularity and Small-world-ness z-scores

The function calculates modularity, number of groups and
small-world-ness z-scores and 99\\ using as null model the list of
networks in the nullDist parameters. Modularity is calculated using the
[`igraph::cluster_spinglass()`](https://r.igraph.org/reference/cluster_spinglass.html)
if the parameter weights is NULL the atribute "weigths" is used, or if
it has the name of an network attribute, that is used as a weigth to
build the modules, when this parameter is NA then no weigth is used.
Only works for one component networks

## Usage

``` r
calc_modularity_swness_zscore(
  g,
  nullDist,
  sLevel = 0.01,
  ncores = 0,
  weights = NA
)
```

## Arguments

- g:

  igraph object

- nullDist:

  list of igraph object with the null model simulations

- sLevel:

  significance level to calculate CI (two tails)

- ncores:

  number of cores to use paralell computation, if 0 sequential
  processing is used.

- weights:

  The weights of the edges. Either a numeric vector or NULL or NA. If it
  is null and the input graph has a ‘weight’ edge attribute then that
  will be used. If NULL and no such attribute is present then the edges
  will have equal weights. Set this to NA if the graph was a ‘weight’
  edge attribute, but you don't want to use it for community detection.

## Value

a list with two data frames: one with indices z-scores and CI

- Clustering:

  Clustering coefficient, measures the average fraction of pairs of
  neighbors of a node that are also neighbors of each other

- PathLength:

  Mean of the shortest paths between all pair of vertices

- Modularity:

  modularity measures how separated are different groups from each
  other, the algorithm `cluster_spinglass` was used to obtain the groups

- zCC,zCP,zMO:

  Z-scores of Clustering,PathLength and Modularity with respect to a
  random Erdos-Renyi null model

- CClow,CChigh,CPlow,CPhigh,MOlow,MOhigh:

  sLevel confidence intervals

- SWness,SWnessCI:

  Small-world-ness and it CI value

- isSW,isSWness:

  Logical variable signalling if the network is Small-world by the
  method of Marina 2018 or the method of Humprhies & Gurney 2008

Another data.frame with the values calculated for the nullDist.

## References

Marina, T. I., Saravia, L. A., Cordone, G., Salinas, V., Doyle, S. R., &
Momo, F. R. (2018). Architecture of marine food webs: To be or not be a
‘small-world.’ PLoS ONE, 13(5), 1–13.
https://doi.org/10.1371/journal.pone.0198217

## Examples

``` r
if (FALSE) { # \dontrun{
nullg <- generateERbasal(netData[[1]],10)
calcModularitySWnessZScore(netData[[1]],nullg)
} # }
```
