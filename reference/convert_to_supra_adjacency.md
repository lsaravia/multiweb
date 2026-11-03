# Convert Multilayer Network to Supra-Adjacency Matrix (SAM)

Convert Multilayer Network to Supra-Adjacency Matrix (SAM)

## Usage

``` r
convert_to_supra_adjacency(
  igraph_list,
  layer_names = NULL,
  interlayer_weight = 1,
  isDirected = TRUE,
  sparse = TRUE,
  use_names = FALSE,
  interlayer = TRUE,
  clean = TRUE
)
```

## Arguments

- igraph_list:

  A list of `igraph` objects (each representing a layer).

- layer_names:

  (Optional) A vector of names for each layer.

- interlayer_weight:

  Numeric, the weight for interlayer links (default: 1).

- isDirected:

  Logical, whether the network is directed.

- sparse:

  Logical, whether to return the supra-adjacency matrix as a sparse
  matrix (default: TRUE).

- use_names:

  Logical; if `TRUE`, edges will use node names instead of numeric IDs
  (default: `FALSE`).

- interlayer:

  Logical; if `FALSE`, no interlayer links are added (default: `TRUE`).

- clean:

  Logical; if `TRUE`, removes empty rows/columns (state nodes with no
  edges) from the output (default: TRUE).

## Value

A list with:

- supra_matrix:

  Supra-adjacency matrix with optional inter-layer links.

- state_nodes_map:

  Mapping of state nodes (layer-node combinations).
