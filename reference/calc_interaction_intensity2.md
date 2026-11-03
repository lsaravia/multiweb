# Calculates the interaction intensity of a food web using allometric scaling relationships to approximate basal metabolic rate

The function uses the body mass of predator/consumer and the number of
prey/resources as source data, then the interaction intensity is
estimated based on O'Gorman (2010) as `a_ij = -b*(M_jˆ0.25 / s_j)`,
where `M_j` is the body mass of predator j, `s_j` is the number of
preys, `b` is a free parameter assumed 0.01. This value of the
interaction strength quantifies the effect of the predator on the prey
per unit of biomass. Assuming a Lotka-Volterra model is equivalent to
the entry A(i,j) of the community matrix, where i is the prey and j the
predator. To predict the direct effect of each prey on its predator,
`a_ji`, we could assume an ecological efficiency, e=0.1, reflecting a
10% transfer of energy between trophic levels, hence `a_ji = e * a_ij`

## Usage

``` r
calc_interaction_intensity2(
  edge_list,
  consumer_n,
  resource_n,
  bodymass_n,
  b = 0.01,
  e = 0.1,
  output_format = "edgelist"
)
```

## Arguments

- edge_list:

  A tibble containing three columns: consumer, resource, and body mass.

- consumer_n:

  The column name for consumers (predators).

- resource_n:

  The column name for resources (prey).

- bodymass_n:

  The column name for body mass of the consumer.

- b:

  Free parameter (default `0.01`).

- output_format:

  Output type: `"edgelist"` (default, containing only predator effects)
  or `"matrix"` (full adjacency matrix).

## Value

A tibble (if `output_format = "edgelist"`) with interaction strengths
**only for predators on prey** on a column named `qRC`. If
`output_format = "matrix"`, an adjacency matrix with interaction
strengths. If `output_format = "igraph" `, an igraph object where edge
weights represent predator effects on prey. Always the IS is a positive
number.

## Details

The function provides two output formats: an edge list with only the
**effect of predators on prey** or an adjacency matrix as a positive
number.

## References

O’Gorman, E. J., Jacob, U., Jonsson, T., & Emmerson, M. C. (2010).
Interaction strength, food web topology and the relative importance of
species in food webs. Journal of Animal Ecology, 79(3), 682–692.
https://doi.org/10.1111/j.1365-2656.2009.01658.x

## Examples

``` r
library(tibble)
#> 
#> Attaching package: ‘tibble’
#> The following object is masked from ‘package:igraph’:
#> 
#>     as_data_frame

# Define an edge list (consumer, resource, bodymass)
edge_list <- tibble(
  predator = c("A", "A", "B", "C"),
  prey = c("B", "C", "D", "D"),
  predator_mass = c(10, 10, 5, 3)  # Body mass for consumers
)

# Compute interaction strength as an edge list (predator effects only)
interaction_strength_edgelist <- calc_interaction_intensity2(
  edge_list, consumer_n = predator, resource_n = prey, bodymass_n = predator_mass, output_format = "edgelist"
)
print(interaction_strength_edgelist)
#> # A tibble: 4 × 3
#>   predator prey      qRC
#>   <chr>    <chr>   <dbl>
#> 1 A        B     0.00281
#> 2 A        C     0.00281
#> 3 B        D     0.00669
#> 4 C        D     0.00760

# Compute interaction strength as an adjacency matrix
interaction_strength_matrix <- calc_interaction_intensity2(
  edge_list, predator, prey, predator_mass, output_format = "matrix"
)
print(interaction_strength_matrix)
#>             A           B           C D
#> A 0.000000000 0.000000000 0.000000000 0
#> B 0.002811707 0.000000000 0.000000000 0
#> C 0.002811707 0.000000000 0.000000000 0
#> D 0.000000000 0.006687403 0.007598357 0
```
