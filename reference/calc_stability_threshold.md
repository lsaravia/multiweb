# Estimate the dynamic stability threshold of an ecological network

This function estimates a stability threshold for an ecological network
using maximum eigenvalue simulations derived from the
[`calc_QSS()`](https://lsaravia.github.io/multiweb/reference/calc_QSS.md)
function. The goal is to approximate the dynamical stability of a system
*without requiring knowledge of the actual interaction strengths or the
real Jacobian matrix*.

## Usage

``` r
calc_stability_threshold(
  g,
  nsim = 1000,
  target_prob = 0.05,
  sd_start = -1,
  sd_step = -0.1,
  sd_min = -20,
  ncores = 0,
  verbose = TRUE
)
```

## Arguments

- g:

  An `igraph` ecological network.

- nsim:

  Integer. Number of QSS simulations per self-damping value.

- target_prob:

  Numeric. Target proportion of negative eigenvalues (default = 0.05).

- sd_start:

  Numeric. Starting value of `selfDamping` (default = -1).

- sd_step:

  Numeric. Increment used to decrease `selfDamping` (default = -0.1).

- sd_min:

  Numeric. Minimum (most negative) value allowed during search.

- ncores:

  Logical. If \> 0 run QSS simulations in parallel.

- verbose:

  Logical. Print diagnostic messages.

## Value

A list with:

- threshold:

  Interpolated `selfDamping` value where the target proportion is
  reached.

- sd_lower:

  Lower bracketing self-damping value.

- sd_upper:

  Upper bracketing self-damping value.

- prob_lower:

  Proportion of negative eigenvalues at `sd_lower`.

- prob_upper:

  Proportion of negative eigenvalues at `sd_upper`.

- sd_values:

  All tested self-damping values.

- probs:

  Corresponding proportions of negative eigenvalues.

- qss_interp:

  Distribution of maximum eigenvalues at the interpolated threshold
  value.

## Details

### Conceptual Background

Real ecosystems persist over ecological timescales, which implies that
their underlying community matrices (Jacobians) must be locally
stable—that is, they produce a negative real part of the maximum
eigenvalue. Because interaction strengths and self-regulation parameters
are rarely known, stability must be inferred indirectly.

The
[`calc_QSS()`](https://lsaravia.github.io/multiweb/reference/calc_QSS.md)
function provides such an indirect method by generating random Jacobian
matrices whose *sign structure* is determined by the interaction
network. Interaction magnitudes are drawn from a uniform distribution
whose limits are defined by the parameters `negative`, `positive`, and
`selfDamping` (Assumed with default values in this function).

The QSS value is the proportion of simulated matrices that are locally
stable (i.e., have negative maximum eigenvalues). A higher QSS implies
that the network topology is more compatible with stable dynamics.

### Purpose of this Function

Because the true amount of self-regulation in a real system is unknown,
this function searches for the value of `selfDamping` for which a target
proportion of negative eigenvalues (typically 5%) is produced. This
"stability threshold" represents the minimal level of self-regulation
required for the network to behave as a persistently stable system.

This approach follows the idea that:

- Real food webs must be dynamically stable.

- Stability requires some level of intraspecific self-regulation.

- The amount of self-regulation needed *depends on network structure*.

By interpolating the crossing point where the QSS proportion reaches a
target value (e.g., 0.05), we obtain an estimate of stability that can
be compared across networks—in the same spirit as structural stability,
modularity–stability relationships, and feasibility analysis.

For details, see References below.

## References

Allesina, S. & Pascual, M. (2008). Network structure, predator–prey
modules, and stability in large food webs. *Theoretical Ecology*, **1**,
55–64.

Grilli, J., Rogers, T. & Allesina, S. (2016). Modularity and stability
in ecological communities. *Nature Communications*, **7**, 12031.

Monteiro, A. B. & Del Bianco Faria, L. (2017). Causal relationships
between population stability and food-web topology. *Functional
Ecology*, **31**, 1294–1300.

Borrelli, J. J. (2015). Selection against instability: stable subgraphs
are most frequent in empirical food webs. *Oikos*, **124**, 1583–1588.

## Examples

``` r
if (FALSE) { # \dontrun{
g <- generate_niche(40, 0.1)
result <- calc_stability_threshold(g, nsim = 500)
plot_stability_curve(result)

# calculate the mean of the negative max eigenvalues at the threshold
mean(result$qss_raw_interp[result$qss_raw_interp < 0])

} # }
```
