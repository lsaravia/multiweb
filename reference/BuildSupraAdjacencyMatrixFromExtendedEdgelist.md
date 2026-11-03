# Build supra-adjacency matrix from edge lists

Build supra-adjacency matrix from edge lists

## Usage

``` r
BuildSupraAdjacencyMatrixFromExtendedEdgelist(
  mEdges,
  Layers,
  Nodes,
  isDirected
)
```

## Arguments

- mEdges:

  data frame with extended edge list, i.e a data frame with columns:
  `node.from, layer.from, node.to, layer.to weight`

- Layers:

  scalar, number of layers

- Nodes:

  scalar, number of nodes

- isDirected:

  logical

## Value

Supra-adjacency matrix, a square matrix of dimension `Nodes * Layers`.
