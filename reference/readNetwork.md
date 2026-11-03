# Read ecological networks in CSV or tab separated file format as edge list or adyacency matrix

If the network is in edge list format and there is a third column is
treated as an edge attribute with the name of the column

## Usage

``` r
readNetwork(
  fileName,
  filePath = NULL,
  fhead = TRUE,
  skipColumn = 1,
  edgeListFormat = 1
)
```

## Arguments

- fileName:

  vector of fileNames with the networks

- filePath:

  path of the files NULL by default

- fhead:

  TRUE if the files have header fields, FALSE otherwise.

- edgeListFormat:

  integer, for the edge list format, if 1 the first column is the "in"
  (the predator), if 2 the first column is the "out" link (the prey).

- skipColum:

  integer, number of columns that are skiped 1 by default

## Value

an igraph object if there is only one file or a list of igraph objects
named after the list without extension

## Examples

``` r

# Reads a network in edge list (interaction list) format, with predators as the first column by default
#
fileName <- system.file("extdata", "WeddellSea_FW.csv", package = "multiweb")
g <- readNetwork(fileName)

# Reads a network in adyacency matrix format, with predators as columns
#
fileName <- system.file("extdata", "BarentsBoreal_FW.csv", package = "multiweb")
g <- readNetwork(fileName)

# Read a vector of files
#
if (FALSE) { # \dontrun{
dn <- list.files("inst/extdata",pattern = "^.*\\.csv$")
netData <- readNetwork(dn,"inst/extdata")
} # }
```
