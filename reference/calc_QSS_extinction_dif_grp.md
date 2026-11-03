# Calculates the QSS difference between the full network and the network minus a group of species

The Quasi-sign stability is estimated with the maximum eigenvalue of the
community matrix (Jacobian) and characterizes the linear stability of
the network. This uses the function
[`calc_QSS()`](https://lsaravia.github.io/multiweb/reference/calc_QSS.md)
so it can take into account the interaction strength if weights are
present. The comparison is made using the Anderson-Darling test with the
function
[`kSamples::ad.test()`](https://rdrr.io/pkg/kSamples/man/ad.test.html)
and the Kolmogorov-Smirnov test
[`stats::ks.test()`](https://rdrr.io/r/stats/ks.test.html), both the
p-values are reported as a measure of strength of the difference.

## Usage

``` r
calc_QSS_extinction_dif_grp(
  g,
  sp_list,
  nsim = 1000,
  ncores = 4,
  istrength = FALSE
)
```

## Arguments

- g:

  igraph network

- sp_list:

  list with the group pf species/nodes we will delete for the comparison

- nsim:

  number of simulations to calculate QSS

- ncores:

  number of cores used to perform the operation

- istrength:

  if TRUE takes the weight attribute of the network as interaction
  strength to calculate QSS.

## Value

a data.frame with:

- Size of the network with deleted nodes

- Connectance of the network with deleted nodes

- Number of components of the network with deleted nodes

  - QSS of the complete network

  - QSS of the network with the deleted nodes

  - difference between the two previous QSS

## Examples

``` r
if (FALSE) { # \dontrun{
g <- netData[[1]]

# Generate random weights
#
V(g)$weight <-  runif(vcount(g))

# Without interaction strength
#
calc_QSS_extinction_dif_grp(g,V(g)$name[1:3],nsim=10,istrength = FALSE)

# With interaction strength
#
calc_QSS_extinction_dif_grp(g,V(g)$name[1:3],nsim=10,istrength = TRUE)
} # }
```
