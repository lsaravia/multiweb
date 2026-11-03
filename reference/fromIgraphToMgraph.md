# From multiple igraph objects to 'mgraph'

This functions takes a list of igraph objects and convert it to a
'mgraph' object.

## Usage

``` r
fromIgraphToMgraph(g, types)
```

## Arguments

- types:

  vector of types that represent the layers with the same length as mg

- mg:

  list of 'igraph' objects

## Value

class 'mgraph' object

## See also

[`readMultiplex()`](https://lsaravia.github.io/multiweb/reference/readMultiplex.md)

## Examples

``` r

# Read a vector of files
#
if (FALSE) { # \dontrun{

fileName <- c(system.file("extdata",  package = "multiweb"))
dn <- list.files("inst/extdata",pattern = "^Kefi2015.*\\.txt$")
g <- readNetwork(dn,"inst/extdata", skipColumn = 2)
fromIgraphToMgraph(g,c("Negative","Positive","Antagonistic"))

} # }
```
