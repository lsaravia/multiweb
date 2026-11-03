# Function to calculate weighted network indices

This function calculates weighted topological indices for ecological
networks, including: weighted linkage density, connectance, generality
and vulnerability and the SDs of the last two based on \[1\] and the
level of omnivory, mean and maximum trophic level based on \[2\] using
[`NetIndices::TrophInd()`](https://rdrr.io/pkg/NetIndices/man/TrophInd.html)
function. The igraph networks must have the weight attribute.

## Usage

``` r
calc_weighted_topological_indices(ig, ncores = 0)
```

## Arguments

- ig:

  An igraph object with the weight attribute representing fluxes

- ncores:

  number of cores used to compute in parallel, if 0 sequential
  processing is used.

## Value

a data.frame with the following fields:

- LD: linkage density

- Connectance: directed Connectance

- TLmean: mean trophic level

- TLmax: maximum trophic level

- LOmnivory: Level of omnivory, quantiﬁes mean of the variety in trophic
  levels of the preys of a consumer

- Vulnerability: mean of number of consumers per prey

- VulSD: the standard deviation of normalized Vulnerability

- Generality: mean number of prey per consumer

- GenSD: the standard deviation of normalized Generality

## References

1.  Bersier, LF. et al. (2002). Quantitative Descriptors Of Food-web
    Matrices. Ecology, 83(9), 2394–2407.

2.  Kones, JK. et al. (2009). Are network indices robust indicators of
    food web functioning? A Monte Carlo approach. Ecological Modelling,
    220, 370–382.

## Examples

``` r

# Generate a test network

g <-   graph_from_literal( 1 -+ 4 -+ 7,2 -+ 5 -+7, 3-+6-+7, 7-+7, 4+-3, simplify = FALSE)

# Add weight

E(g)$weight <- sample(c(.1,.2,.8,.9),gsize(g),replace=TRUE)

calc_weighted_topological_indices(g)
#>        wLD wConnectance  wTLmean   wTLmax wLOmnivory wVulnerability    wVulSD
#> 1 1.839137    0.2627338 1.626374 2.384615 0.08140971       1.083482 0.5894277
#>   wGenerality   wGenSD
#> 1    2.594791 2.067256
```
