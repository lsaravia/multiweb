# Calculate the incoherence index of a food web

The incoherence index is based in how the species fit in discrete
trophic levels when Q is closer to 0 more coherent and stable is a food
web, and the less omnivory it has. It calculates the trophic level using
the package NetIndices.

## Usage

``` r
calc_incoherence(ig, ncores = 0)
```

## Arguments

- ig:

  an igraph object or a list of igraph objects

- ncores:

  number of cores used to compute in parallel, if 0 sequential
  processing is used.

## Value

a data.frame with the following fields

- Q:

  incoherence (0=coherent)

- rQ:

  ratio of Q with expected Q under null expectation of a random network
  given N=nodes L=links B=basal nodes Lb=basal links

- mTI:

  mean trophic level

- rTI:

  ratio of mTI with expected TI under the same null model expectation
  than Q

## Details

Based on:

## References

Johnson, S., Domínguez-García, V., Donetti, L., & Muñoz, M. A. (2014).
Trophic coherence determines food-web stability. Proceedings of the
National Academy of Sciences , 111(50), 17923–17928.
https://doi.org/10.1073/pnas.1409077111

Johnson, S., & Jones, N. S. (2017). Looplessness in networks is linked
to trophic coherence. Proceedings of the National Academy of Sciences,
114(22), 5618–5623. https://doi.org/10.1073/pnas.1613786114

## Examples

``` r

calc_incoherence(netData[[1]])
#> Warning: executing %dopar% sequentially: no parallel backend registered
#> 
#> Attaching package: ‘igraph’
#> The following objects are masked from ‘package:stats’:
#> 
#>     decompose, spectrum
#> The following object is masked from ‘package:base’:
#> 
#>     union
#> Loading required package: MASS
#>           Q        rQ      mTI       rTI
#> 1 0.5248904 0.2063463 2.897171 0.3650221

```
