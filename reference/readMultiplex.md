# Read ecological multiplex networks using different files for each layer as a 'node-aligned' multilayer network

This functions uses
[`readNetwork()`](https://lsaravia.github.io/multiweb/reference/readNetwork.md)
to read files that represent network layers that can be different
interaction types, represented by the type attribute of the igraph
object.

## Usage

``` r
readMultiplex(
  fileName,
  types = c("Negative", "Positive", "Antagonistic"),
  filePath = NULL,
  fhead = TRUE,
  skipColumn = 1,
  format = "layers"
)
```

## Arguments

- fileName:

  vector of fileNames with the layers of the networks

- types:

  vector of types that represent the layers

- filePath:

  path of the files NULL by default

- fhead:

  TRUE if the files have header fields, FALSE otherwise.

- format:

  string, "layers" is the default were different layers are coded as
  different files there must be the same number of files as the length
  of the types vector as each type represent a layer. "GLV" Represent
  multiple interaction types as pairs of entries in a matrix so
  competition is represented as `a[i,j] = -1, a[j,i] = -1`, predation
  `a[i,j] = 1, a[j,i] = -1` where species j is the predator. Mutualism
  is `a[i,j] = a[j,i] = 1`. Any negative or positive number works
  because intensity of interactions are not registered.

- skipColum:

  integer, number of columns that are skiped 1 by default

## Value

an mgraph object that is a list of igraph objects with the type
attribute set to each kind of layer or interaction

## See also

[`readNetwork()`](https://lsaravia.github.io/multiweb/reference/readNetwork.md)

## Examples

``` r

# Read a vector of files
#
if (FALSE) { # \dontrun{

fpath <- system.file("extdata", package = "multiweb")
dn <- list.files(fpath,pattern = "^Kefi2015.*\\.txt$")
netData <- readMultiplex(dn,c("Negative","Positive","Antagonistic"),fpath,skipColum=2)
} # }
```
