# Calculation Small-world-ness z-scores

The function calculates small-world-ness z-scores and 99\\ using as null
model the list of networks in the nullDist parameter.

## Usage

``` r
calc_swness_zscore(g, nullDist, sLevel = 0.01, ncores = 0, weights = NA)
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

  Z-scores of Clustering,PathLength and Modularity with respect to
  `nullDist` models

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
calc_swness_zscore(netData[[1]],nullg)
} # }
```
