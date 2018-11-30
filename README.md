# EcoNetwork 

This is an R package to calculate recently developed network metrics and also some old ones conveniently in the same place. Functions to read, plot and a big set of networks are also included in the package. 

## Instalation 

```R
require(devtools)
install_github("lsaravia/EcoNetworks")
```

## Examples

Read the example data included in the package

```R
require(EcoNetwork)

fileName <- c(system.file("extdata",  package = "EcoNetwork"))
fileName <- system.file("extdata", "BarentsBoreal_FW.csv", package = "EcoNetwork")
g <- readNetwork(fileName)


```



## References

1. Marina, T. I., Saravia, L. A., Cordone, G., Salinas, V., Doyle, S. R., & Momo, F. R. (2018). Architecture of marine food webs: To be or not be a ‘small-world.’ PLoS ONE, 13(5), 1–13. https://doi.org/10.1371/journal.pone.0198217

2. Saravia, L. A., Marina, T. I., De Troch, M., & Momo, F. R. (2018). Ecological Network assembly: how the regional meta web influence local food webs. BioRxiv, 340430. Retrieved from http://biorxiv.org/content/early/2018/06/21/340430.abstract

