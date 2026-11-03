# PageRank importance of species in a food web following Allesina & Pascual (2009)

Computes the "food-web PageRank" of Allesina & Pascual (2009). Unlike
the standard PageRank applied directly to a food web (e.g.
[`igraph::page_rank()`](https://r.igraph.org/reference/page_rank.html)
with links from prey to predator), which accumulates importance in top
predators, this version measures how much other species *depend* on each
species, so that importance flows from consumers to their resources.

## Usage

``` r
calc_pagerank_ap(g, weights = FALSE, remove_loops = TRUE)
```

## Arguments

- g:

  igraph directed food web with links from prey to predator (the
  multiweb convention).

- weights:

  logical, if TRUE and `g` has a `weight` edge attribute, each consumer
  splits its importance proportionally to the link weights instead of
  equally. The links added to/from the root have weight 1 for basal
  species and, for the detritus links, the mean weight of the network.
  The original method is unweighted (default FALSE).

- remove_loops:

  logical, remove cannibalistic self-loops (default TRUE).

## Value

A data frame sorted by decreasing importance with columns:

- Species:

  species name (or vertex index if the graph has no names).

- PageRank:

  importance, summing 1 over species (root excluded).

- Rank:

  rank of the species (1 = most important; ties get the minimum rank).

- Basal:

  TRUE for species without prey, that are linked to the root.

## Details

The algorithm:

1.  Adds a "root" node (external energy/nutrient source) that is the
    only resource of the basal species (species without prey).

2.  Adds a link from every species to the root, representing loss of
    matter (detritus/recycling). This makes the graph irreducible and
    primitive, so no damping factor is needed (damping = 1), as argued
    by the authors.

3.  Each consumer splits its importance equally (or proportionally to
    link weights, if `weights = TRUE`) among its resources; the
    importance vector is the dominant eigenvector of the resulting
    column-stochastic matrix.

The root is removed from the output and the scores of the species are
rescaled to sum 1.

The method is defined for trophic networks only: the root node
represents the flow of matter, which has no meaning for non-trophic
layers.

## References

Allesina, S. & Pascual, M. (2009). Googling food webs: can an
eigenvector measure species' importance for coextinctions? PLoS
Computational Biology, 5(9), e1000494.
[doi:10.1371/journal.pcbi.1000494](https://doi.org/10.1371/journal.pcbi.1000494)

## See also

[`calc_centrality()`](https://lsaravia.github.io/multiweb/reference/calc_centrality.md)

## Examples

``` r
# Linear chain: resource -> herbivore -> predator
g <- igraph::graph_from_literal(Resource -+ Herbivore, Herbivore -+ Predator)
calc_pagerank_ap(g)
#>     Species  PageRank Rank Basal
#> 1  Resource 0.5000000    1  TRUE
#> 2 Herbivore 0.3333333    2 FALSE
#> 3  Predator 0.1666667    3 FALSE

# Compare with the standard PageRank, which ranks the predator first
igraph::page_rank(g)$vector
#>  Resource Herbivore  Predator 
#> 0.1844168 0.3411710 0.4744122 
```
