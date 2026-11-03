# Calculates the quantitative connectance, and the effective number of flows using the interaction matrix and a vector of species abundances

The Quantitative connectance (Cq) takes into account both the
distribution of per capita interaction strengths among species in the
web and the distribution of species’ abundances and quantifies the
diversity of network fluxes, if all the species have the same flux is
equal to the directed connectance. The mean or effective number of flows
impinging upon or emanating from a tipical node (LDq) is based the
average flow diversity, and when all flows are equal is similar to
linkage density. Both measures are based in Shannon information theory.
The total interaction flux is measured as
`T[i,j] <- d[i] * d[j] * interM[i,j]`. The effective Cq is calculated
following the formulas in appendix 2 of \[1\], LDq follows \[2\]

## Usage

``` r
calc_quantitative_connectance(interM, d)
```

## Arguments

- interM:

  per capita interaction strength matrix

- d:

  species' abundances vector

## Value

A list with Cq,the quantitative connectance index and LDq, =

## References

1.  Fahimipour, A.K. & Hein, A.M. (2014). The dynamics of assembling
    food webs. Ecol. Lett., 17, 606–613

2.  Ulanowicz, R.E. & Wolff, W.F. (1991). Ecosystem flow networks:
    Loaded dice? Math. Biosci., 103, 45–68

## Examples

``` r

# 3 predators 2 preys unequal fluxes
#
m <- matrix()
matrix(0,nrow=4,ncol=4)
#>      [,1] [,2] [,3] [,4]
#> [1,]    0    0    0    0
#> [2,]    0    0    0    0
#> [3,]    0    0    0    0
#> [4,]    0    0    0    0
m[1,2] <- m[1,3] <- m[3,4]<- .2
#> Error in `[<-`(`*tmp*`, 3, 4, value = 0.2): subscript out of bounds
m[2,1] <- m[3,1] <- m[4,3] <- -2
#> Error in `[<-`(`*tmp*`, 4, 3, value = -2): subscript out of bounds

calc_quantitative_connectance(m, c(1,1,1,1))
#> Error in calc_quantitative_connectance(m, c(1, 1, 1, 1)): length of species vector has to be equal to interM dimensions

# Equal input and output fluxes

m <- matrix(0,nrow=4,ncol=4)
m[1,2] <- m[1,3] <- m[3,4]<- 2
m[2,1] <- m[3,1] <- m[4,3] <- -2
calc_quantitative_connectance(m, c(1,1,1,1))
#> $Cq
#> [1] 0.375
#> 
#> $LDq
#> [1] 1.587401
#> 
```
