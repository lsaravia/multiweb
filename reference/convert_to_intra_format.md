# Convert a List of igraph Objects to Intralayer Edge Format

This function converts a multilayer network, represented as a list of
igraph objects, into a standardized format suitable for Infomap and
other multilayer network analyses.

## Usage

``` r
convert_to_intra_format(igraph_list, use_names = FALSE)
```

## Arguments

- igraph_list:

  A list of `igraph` objects, each representing a network layer.

- use_names:

  Logical; if `TRUE`, edges will use node names instead of numeric IDs
  (default: `FALSE`).

## Value

A list with two data frames:

- vertices:

  A data frame containing node_id and corresponding node names.

- intra:

  A data frame containing intra-layer edges with columns: `layer_id`,
  `node_id1`, `node_id2`, `weight`.

## Details

Each node is assigned a unique numeric ID, and edges across layers are
recorded with layer-specific identifiers.

## Examples

``` r
if (FALSE) { # \dontrun{
  multi_format <- convert_to_intra_format(list(layer1, layer2))
  print(multi_format$vertices)
  print(multi_format$intra)
} # }
```
