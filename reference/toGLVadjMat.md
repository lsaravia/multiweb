# From multiple interaction object 'mgraph'to GLV adjacency matrix

This functions takes a 'mgraph' object and convert it to a Generalized
Lotka-Volterra adjacency matrix Position is important in fact the order
is Competitive/Negative,Mutualistic/Positive,Trophic/Antagonistic If
'istrength' is TRUE the attribute weight is assumed as the strength of
the interacion

## Usage

``` r
toGLVadjMat(
  mg,
  types = c("Negative", "Positive", "Antagonistic"),
  istrength = FALSE
)
```

## Arguments

- mg:

  multiple interaction object, class 'mgraph'

- types:

  vector of types that represent the layers

- istrength:

  edge weights

## Value

an matrix with the signs of the interactions

## See also

[`readMultiplex()`](https://lsaravia.github.io/multiweb/reference/readMultiplex.md)

## Examples

``` r

# Read a vector of files
#
if (FALSE) { # \dontrun{

fpath <- system.file("extdata", package = "multiweb")
dn <- list.files(fpath,pattern = "^Kefi2015.*\\.txt$")
netData <- readMultiplex(dn,c("Negative","Positive","Antagonistic"),fpath,skipColum=2)
toGLVadjMat(netData)
} # }
```
