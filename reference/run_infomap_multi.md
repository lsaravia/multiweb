# Run Infomap on Multilayer Networks with Intralayer Links

This function exports a multilayer network (as a list of igraph objects)
into Infomap's required format, runs Infomap via an external binary, and
imports the detected communities back into R as a data frame. It
requires Infomap to be installed on the system. See installation
instructions at [Infomap](https://www.mapequation.org/infomap/#Install).

## Usage

``` r
run_infomap_multi(
  igraph_list,
  layer_names = NULL,
  infomap_path = "infomap",
  output_dir = tempdir(),
  directed = TRUE,
  two_level = TRUE,
  multilayer_relax_rate = 0.15,
  seed = 123
)
```

## Arguments

- igraph_list:

  A list of `igraph` objects, each representing a network layer.

- layer_names:

  A character vector with layer names corresponding to each igraph
  object (default: `NULL`).

- infomap_path:

  Character string specifying the path to the Infomap binary (default:
  `"infomap"`, assumes it is in the system `PATH`).

- output_dir:

  Character string specifying the directory for Infomap results
  (default: [`tempdir()`](https://rdrr.io/r/base/tempfile.html)).

- directed:

  Logical; if `TRUE`, the network is treated as directed (default:
  `TRUE`).

- two_level:

  Logical; if `TRUE`, uses two-level Infomap instead of hierarchical
  Infomap (default: `TRUE`).

- multilayer_relax_rate:

  Numeric; relaxation rate for multilayer links (default: `0.15`).

- seed:

  Numeric; random seed for Infomap to ensure reproducibility (default:
  `123`).

## Value

A data frame containing the detected modules with columns:

- module:

  Module/community assignment from Infomap.

- node:

  The node name from the original igraph objects.

- layer:

  The corresponding layer name from `layer_names`.

- flow:

  The fraction of flow assigned to the module.

## Details

If the network has the `weight` attribute, it will be used as the weight
of the edges. Since the exported multilayer network structure assumes
only `*Intra` links (i.e., no explicit inter-layer links), Infomap
internally constructs inter-layer connections for each physical node
using the `--multilayer-relax-rate` parameter (default: 0.15). This
creates inter-layer links between all instances of the same node across
layers. The weights of these links are proportional to the weighted
out-degree of the node in the target layer, resulting in non-uniform
inter-layer transition probabilities.

This implicit construction of inter-layer links enables Infomap to model
node-aligned multiplex dynamics without requiring an explicit `*Inter`
section in the input file.

## References

D. Edler, A. Holmgren and M. Rosvall, *The MapEquation software
package*, available online at <https://www.mapequation.org>.

## Examples

``` r
if (FALSE) { # \dontrun{
# Load network data
fileName <- system.file("extdata", package = "multiweb")
dn <- list.files(fileName, pattern = "^Kefi2015.*\\.txt$")
g <- readNetwork(dn, fileName, skipColumn = 2)
class(g)
names(g) <- c("Negative", "Positive", "Trophic")

# Run Infomap
run_infomap_multi(g)
} # }
```
