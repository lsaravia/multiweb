# Calculate Topological Roles Within Network Modules

This function identifies the topological roles of nodes (e.g., species)
in an ecological or interaction network based on their position within
and among modules. Modules (or communities) are detected using the
[`cluster_spinglass`](https://r.igraph.org/reference/cluster_spinglass.html)
algorithm by default, unless an alternative community object is
provided.

## Usage

``` r
calc_topological_roles(g, nsim = 1, ncores = 0, community = NULL)
```

## Arguments

- g:

  An `igraph` object representing the network.

- nsim:

  Integer. Number of independent runs of the Spinglass community
  detection algorithm (default is 1).

- ncores:

  Integer. Number of cores to use for parallel computation. Set to 0 for
  sequential processing.

- community:

  Optional. A community detection object (class `communities`) from an
  alternative algorithm. If `NULL`,
  [`cluster_spinglass`](https://r.igraph.org/reference/cluster_spinglass.html)
  is used.

## Value

A data frame with the following columns:

- node:

  Node ID.

- within_module_degree:

  Within-module degree z-score (\\z\\).

- among_module_conn:

  Participation coefficient (\\PC\\).

## Details

Topological roles are quantified using two parameters:

- \\z\\ (within-module degree z-score): quantifies how well-connected a
  node is to others within its own module, relative to other members of
  the same module.

- \\PC\\ (participation coefficient): quantifies how evenly a node’s
  links are distributed across different modules.

When using the Spinglass algorithm, which is stochastic, multiple runs
can be performed by setting the `nsim` parameter to repeat the
procedure.

## References

1.  Guimerà, R. & Nunes Amaral, L.A. (2005). Functional cartography of
    complex metabolic networks. *Nature*, 433, 895–900.

2.  Kortsch, S., Primicerio, R., Fossheim, M., Dolgov, A.V., &
    Aschan, M. (2015). Climate change alters the structure of Arctic
    marine food webs due to poleward shifts of boreal generalists.
    *Proceedings of the Royal Society B: Biological Sciences*,
    282(1814), 20151546.
    [doi:10.1098/rspb.2015.1546](https://doi.org/10.1098/rspb.2015.1546)

## Examples
