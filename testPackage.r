# Test

fileName <- c(system.file("extdata",  package = "EcoNetwork"))
fileName <- system.file("extdata", "BarentsBoreal_FW.csv", package = "EcoNetwork")
g <- readEcoNetwork(fileName)



dn <- list.files("inst/extdata",pattern = "^.*\\.csv$")
netData <- readEcoNetwork(dn,"inst/extdata")


plotEcoNetworkTrophLevel(netData[[1]])

require(igraph)
g <-   graph_from_literal( 1 -+ 4 -+ 7,2 -+ 5 -+7, 3-+6-+7, 7-+7, 4+-3, simplify = FALSE)
topologicalIndices(g)
plotTrophLevel(g,vertexLabel = TRUE,vertexSizeFactor = 20)
calcIncoherence(g)

g <-   graph_from_literal( 1 -+ 4 -+ 7,2 -+ 5 -+7, 3-+6-+7, 4+-5, simplify = FALSE)
topologicalIndices(g)
plotTrophLevel(g,vertexLabel = TRUE,vertexSizeFactor = 20)
calcIncoherence(g)


