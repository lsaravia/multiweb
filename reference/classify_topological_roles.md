# Classify and plot topological roles

This function classifies species into four topological roles based on
their within-module degree (dz) and among-module connectivity (PC),
following the approach of Kortsch (2015). The classification is as
follows:

- **modhub**: Species with high within-module degree (dz \>= 2.5) and
  low among-module connectivity (PC \< 0.625).

- **hubcon**: Species with high within-module degree (dz \>= 2.5) and
  high among-module connectivity (PC \>= 0.625).

- **modspe**: Species with low within-module degree (dz \< 2.5) and low
  among-module connectivity (PC \< 0.625).

- **modcon**: Species with low within-module degree (dz \< 2.5) and high
  among-module connectivity (PC \>= 0.625).

## Usage

``` r
classify_topological_roles(tRoles, g, community = NULL)
```

## Arguments

- tRoles:

  A data frame of calculated topological roles obtained from
  [`calc_topological_roles()`](https://lsaravia.github.io/multiweb/reference/calc_topological_roles.md).

- g:

  An igraph network object.

- community:

  An igraph community object used to determine module membership for
  plotting. If NULL, the function calculates community structure using
  `cluster_spinglass()`.

## Value

A list containing:

- `hub_conn`: A data frame with node classifications and topological
  role values.

- `gp`: A ggplot object visualizing the topological roles.

## Details

The function also generates a ggplot visualization of the species'
topological roles, with nodes colored according to their module
membership.

## References

1.  Kortsch, S., Primicerio, R., Fossheim, M., Dolgov, A. V., &
    Aschan, M. (2015). Climate change alters the structure of arctic
    marine food webs due to poleward shifts of boreal generalists.
    Proceedings of the Royal Society B: Biological Sciences,
    282(1814), 20151546. https://doi.org/10.1098/rspb.2015.1546

2.  Saravia, L. A., Marina, T. I., Kristensen, N. P., De Troch, M., &
    Momo, F. R. (2022). Ecological network assembly: How the regional
    metaweb influences local food webs. Journal of Animal Ecology,
    91(3), 630–642. https://doi.org/10.1111/1365-2656.13652

## Examples

``` r
if (FALSE) { # \dontrun{

g <- netData[[2]]

tp <- calc_topological_roles(g,nsim=10)

classify_topological_roles(tp,g)
} # }
```
