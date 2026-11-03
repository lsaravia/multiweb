# Harmonize Node Sets Across Layers in a Multiplex Network

This function ensures that all igraph layers in a multiplex network have
the same set of nodes. Nodes missing in any layer are added as isolated
nodes (no edges).

## Usage

``` r
harmonize_node_sets(g.list)
```

## Arguments

- g.list:

  A list of igraph objects representing the multiplex layers.

## Value

A list of igraph objects with identical node sets.

## Examples

``` r
g.list <- harmonize_node_sets(g.list)
#> Error: object 'g.list' not found
```
