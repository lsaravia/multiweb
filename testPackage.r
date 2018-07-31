# Test

fileName <- c(system.file("extdata",  package = "EcoNetwork"))
fileName <- system.file("extdata", "BarentsBoreal_FW.csv", package = "EcoNetwork")
g <- readEcoNetwork(fileName)



dn <- list.files("inst/extdata",pattern = "^.*\\.csv$")
netData <- readEcoNetwork(dn,"inst/extdata")


plotEcoNetworkTrophLevel(netData[[1]])

require(igraph)
require(NetIndices)
g <-   graph_from_literal( 1 -+ 4 -+ 7,2 -+ 5 -+7, 3-+6-+7, 7-+7, 4+-3, simplify = FALSE)
topologicalIndices(g)
plotTrophLevel(g,vertexLabel = TRUE,vertexSizeFactor = 20)
calcIncoherence(g)

g <-   graph_from_literal( 1 -+ 4 -+ 7,2 -+ 5 -+7, 3-+6-+7, 4+-5, simplify = FALSE)
topologicalIndices(g)
plotTrophLevel(g,vertexLabel = TRUE,vertexSizeFactor = 20)
calcIncoherence(g)
calcModularitySWnessZScore(g,10,0.1,paralell = FALSE)

# test reading .dat files
#

dn <- list.files("~/Dropbox/Projects/NetworksAsUnifyingPrinciple/Data",pattern = "^.*\\.dat$")
nets <- readNetwork(dn,"~/Dropbox/Projects/NetworksAsUnifyingPrinciple/Data",fhead=FALSE)
topologicalIndices(list(netData$CaymanIs_FW,nets$cayman_islands))
require(RColorBrewer)
plotTrophLevel(nets$cayman_islands)
plotTrophLevel(netData$CaymanIs_FW)

# The incoherence parameter seems not well calculated compared to Johnson 2017 Appendix Table S1
#
require(igraph)

calcIncoherence(nets$benguela)
mean(degree(nets$benguela))

calcIncoherence(nets$bridge)
mean(degree(nets$benguela))


# Test modularity
#
g <-   graph_from_literal( 1 -+ 4 -+ 7,2 -+ 5 -+ 9, 4+-5,
                           3 -+ 6 -+ 8,5 -+8, simplify = FALSE)
topologicalIndices(g)
plotTrophLevel(g,vertexLabel = TRUE,vertexSizeFactor = 20,modules = TRUE)
title(main="9 species")
calcIncoherence(g)
calcModularitySWnessZScore(g,100,0.1,paralell = FALSE)

g <-   graph_from_literal( 1 -+ 4 -+ 7,2 -+ 5 -+ 9, 4+-5,
                           3 -+ 6 -+ 8,5 -+8, 8-+ 9, simplify = FALSE)
plotTrophLevel(g,vertexLabel = TRUE,vertexSizeFactor = 20,modules = TRUE)
calcIncoherence(g)
calcModularitySWnessZScore(g,100,0.1,paralell = FALSE)
