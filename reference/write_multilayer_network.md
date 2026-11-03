# Write a Multilayer Network to a File in Infomap-Compatible Format

This function writes a multilayer network (represented as a list of
igraph objects) to a text file in Infomap's required \*Intra format.

## Usage

``` r
write_multilayer_network(igraph_list, file_path)
```

## Arguments

- igraph_list:

  A list of `igraph` objects, each representing a layer.

- file_path:

  A string specifying the path to the output file.

## Value

A list with two data frames:

- vertices:

  A data frame of node IDs and names.

- intra:

  A data frame of intra-layer edges with columns: `layer_id`,
  `node_id1`, `node_id2`, `weight`.

## Details

The output file consists of:

- A list of unique vertices with their numeric IDs.

- Intra-layer edges with weights.

## Examples

``` r
if (FALSE) { # \dontrun{
  write_multilayer_network(list(layer1, layer2), "network.net")
} # }
```
