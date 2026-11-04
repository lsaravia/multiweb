# Interface R and igraph with Infomap (external binary)

This function exports an igraph network to Infomap's required format,
runs Infomap, and imports the detected communities back into R as an
igraph-compatible cluster object. It requires the installation of
infomap in the system from https://www.mapequation.org/infomap/#Install
if the network has the weight attribute, it will be used as the weight
of the edges. The modularity stored in the returned object is computed
from the partition using the weights (if present) and the `directed`
argument. It works for networks with more than 1 connected component,
but not for networks with isolated nodes.

## Usage

``` r
run_infomap(
  graph,
  infomap_path = "infomap",
  output_dir = tempdir(),
  directed = TRUE,
  two_level = TRUE,
  seed = 123,
  return_df = FALSE
)
```

## Arguments

- graph:

  An igraph object.

- infomap_path:

  Path to the Infomap binary (default assumes it's in system PATH).

- output_dir:

  Temporary directory for Infomap results.

- directed:

  Boolean: Treat the network as directed (default TRUE).

- two_level:

  Boolean: Use two-level Infomap (default TRUE).

- seed:

  Numeric: Random seed for Infomap (default 123).

- return_df:

  Boolean: If TRUE, return a data.frame with `node`, `module`, and
  `flow`; if FALSE, return a community object (default).

## Value

Either an igraph community object or a list with a data.frame and
codelength.

## Examples

``` r
# Example with the intermediate files in actual folder
py_infomap <- run_infomap(netData[[1]], output_dir = ".")

membership(py_infomap)
#>        Trachurus       Sardinella       Sciaenidae          Ariidae 
#>                3                3                3                1 
#>       Merluccius           Dentex   Pseudotolithus       Serranidae 
#>                3                1                1                3 
#>      Polynemidae  P. longirostris       Haemulidae Aristeu varidens 
#>                3                3                1                1 
#>   Small pelagics     Mesopelagics   Large pelagics    Demersal fish 
#>                3                3                2                3 
#>      Condrichtis     Cephalopodes     Miscelaneous          Penaeus 
#>                2                1                3                3 
#>   Marine mammals     Marine birds     Macrobenthos    Middlebenthos 
#>                4                1                3                3 
#>        Large Zoo        Small Zoo    Phytoplankton         Detritus 
#>                3                3                3                3 
modularity(py_infomap)
#> [1] 0.07638415
```
