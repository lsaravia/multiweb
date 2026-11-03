# Calc the Quasi Sign Stability measure for antagonistic (predator-prey) or mgraph networks with multiple interactions.

The QSS measure is the proportion of matrices that are locally stable,
these matrices are created by sampling the values of the community
matrix (the Jacobian) from a uniform distribution, preserving the sign
structure \[1\]. If the 'ig' parameter is an `mgraph` network it needs
to have been built with the order
`c("Competitive", "Mutualistic", "Trophic")` It also calculates the mean
of the real part of the maximum eigenvalue, which is also a measure of
stability \[2\]. It uses a uniform distribution between 0 and maximum
values given by the parameters `negative`, `positive` and `selfDamping`,
corresponding to the sign of interactions and self-limitation effect
\[3,4\]. If the edges of the networks have a weight attribute and
`istrength` parameter is true, weight will be used as interaction
strength, then the limits of the uniform distribution will be
`negative*-x`, `positive*x`, `selfDamping*x`, where x is the value of
the weigth for the edge. If the values of these parameters are 0 then
there is no interaction of that kind. The default values for `negative`,
`positive` and `selfDumping` assume a maximum ecological transfer
efficiency of 10%.

## Usage

``` r
calc_QSS(
  ig,
  nsim = 1000,
  ncores = 0,
  negative = -10,
  positive = 1,
  selfDamping = -1,
  istrength = FALSE,
  returnRaw = FALSE
)
```

## Arguments

- ig:

  igraph or a list of igraph networks or mgraph network

- nsim:

  number of simulations to calculate QSS, if the number of simulations
  is 1 then it calculates the maximum eingenvalue for the mean of
  interaction strength, if `istrength==FALSE` the mean interaction
  strength are `negative/2`, `positive/2`, `selfDamping/2`.

- ncores:

  number of cores to use in parallel computation if 0 do not set future
  strategy inside function

- negative:

  the maximum magnitude of the negative interaction (the effect of the
  predator on the prey) must be \<= 0

- positive:

  the maximum magnitude of the positive interaction (the effect of the
  prey on the predator) must be \>= 0

- selfDamping:

  the maximum magnitude of the self-limitation (the effect of the
  species on itself) must be \<= 0, if = 0 or \> 0 then the value is a
  random number between 0 and `selfDamping`.

- istrength:

  If TRUE takes the weigth attribute of the network as interaction
  strength.

- returnRaw:

  if TRUE returns all the values of the maximum eingenvalues

## Value

if parameter `returnRaw` is `FALSE` returns a data.frame with the QSS,
and MEing, the mean of the real part of the maximum eingenvalue. If
`returnRaw` is `TRUE` it returns the values of randomized real part of
maximum eingenvalue (maxre)

## References

1.  Allesina, S. & Pascual, M. (2008). Network structure, predator -
    Prey modules, and stability in large food webs. Theor. Ecol., 1,
    55–64.

2.  Grilli, J., Rogers, T. & Allesina, S. (2016). Modularity and
    stability in ecological communities. Nat. Commun., 7, 12031

3.  Monteiro, A.B. & Del Bianco Faria, L. (2017). Causal relationships
    between population stability and food-web topology. Functional
    Ecology, 31, 1294–1300.

4.  Borrelli, J. J. 2015. Selection against instability: stable
    subgraphs are most frequent in empirical food webs. - Oikos 124:
    1583–1588.

## Examples

``` r
if (FALSE) { # \dontrun{

g <- netData[[2]]

tp <- calc_QSS(g)

# Read Multiplex network and calculate QSS

fileName <- c(system.file("extdata",  package = "multiweb"))
dn <- list.files(fileName,pattern = "^Kefi2015.*\\.txt$")
gt <- readMultiplex(dn,types,"inst/extdata", skipColumn = 2)
calc_QSS(gt)

} # }
```
