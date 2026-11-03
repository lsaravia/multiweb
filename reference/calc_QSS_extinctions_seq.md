# Calculates the QSS for an extinction sequence

This functions calculates the QSS
[`calc_QSS()`](https://lsaravia.github.io/multiweb/reference/calc_QSS.md)
for a sequence of incremental deletions of species (nodes) given by the
vector `seq`, it does not produce secondary extinctions

## Usage

``` r
calc_QSS_extinctions_seq(
  g_del,
  seq,
  nsim = 1000,
  ncores = 4,
  istrength = FALSE
)
```

## Arguments

- g_del:

  igraph object for the deletion sequence

- seq:

  vector or list with the nodes for the extinction sequence, the
  extinctions are incremental.

- nsim:

  number of simulations used in QSS calculation

- ncores:

  number of cores used to calculate QSS

- istrength:

  parameter istrength for QSS function

## Value

data frame with the number of remaining nodes, the connectance, the
number unconnected of components, the median QSS, and the name of the
last deleted node.

## Examples

``` r
if (FALSE) { # \dontrun{
g <- netData[[1]]

# Generate random weights
#
V(g)$weight <-  runif(vcount(g))

# Without interaction strength
#
calc_QSS_extinctions_seq(g,V(g)$name[1:3],nsim=10,istrength = FALSE)

# With interaction strength
#
calc_QSS_extinctions_seq(g,V(g)$name[1:3],nsim=10,istrength = TRUE)
} # }
```
