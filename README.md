
# multiweb 

[![DOI](https://zenodo.org/badge/142077196.svg)](https://zenodo.org/badge/latestdoi/142077196)

This is an R package to analize ecological networks, including multilayer networks, recently developed network metrics and also some old ones conveniently in the same place. Functions to read, plot and a big set of networks are also included in the package. 

## Instalation 

```R
require(devtools)
install_github("lsaravia/multiweb")
```

## Examples

Read the example data included in the package

```R
require(multiweb)

fileName <- c(system.file("extdata",  package = "multiweb"))
fileName <- system.file("extdata", "BarentsBoreal_FW.csv", package = "multiweb")
g <- readNetwork(fileName)


```

Read multiple interaction network in different layers

```R

fileName <- c(system.file("extdata",  package = "multiweb"))
dn <- list.files("inst/extdata",pattern = "^Kefi2015.*\\.txt$")
g <- readNetwork(dn,"inst/extdata", skipColumn = 2)
gt <- igraph2mgraph(g,c("Negative","Positive","Antagonistic"))

types <- c("Competitive","Mutualistic","Trophic")
gt <- readMultiplex(dn,types,"inst/extdata", skipColumn = 2)

```

Calculate QSS (quasi-sign-stability) for multiple interactions networks (mgraph)

```R

calc_QSS(gt)


```



## References

1. Marina, T. I., Saravia, L. A., Cordone, G., Salinas, V., Doyle, S. R., & Momo, F. R. (2018). Architecture of marine food webs: To be or not be a ‘small-world.’ PLoS ONE, 13(5), 1–13. https://doi.org/10.1371/journal.pone.0198217

2. Ecological Network assembly: how the regional metaweb influences local food webs
Leonardo A. Saravia, Tomás I. Marina, Marleen De Troch, Fernando R. Momo
bioRxiv 340430; doi: https://doi.org/10.1101/340430
