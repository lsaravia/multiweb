# Test

fileName <- c(system.file("extdata",  package = "EcoNetwork"))
fileName <- system.file("extdata", "BarentsBoreal_FW.csv", package = "EcoNetwork")
g <- readNetwork(fileName)
g


dn <- list.files("inst/extdata",pattern = "^.*\\.csv$")
g<- readNetwork(dn,"inst/extdata")
names(g)
E(g[[1]])[which_loop(g[[1]])]

plotTrophLevel(g[[1]])

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


# test reading .csv files
#
dn <- list.files("~/Academicos/GitProjects/MarineFoodWebsSmallWorld/Data/FoodWebs",pattern = "^.*\\.csv$")
#
# BEWARE different files have different structure
#
nets <- readNetwork(dn[11],"~/Academicos/GitProjects/MarineFoodWebsSmallWorld/Data/FoodWebs")
nets <- readNetwork(dn[16],"~/Academicos/GitProjects/MarineFoodWebsSmallWorld/Data/FoodWebs")
nets <- readNetwork(dn[17],"~/Academicos/GitProjects/MarineFoodWebsSmallWorld/Data/FoodWebs")
nets <- readNetwork(dn[18],"~/Academicos/GitProjects/MarineFoodWebsSmallWorld/Data/FoodWebs")
nets <- readNetwork(dn[19],"~/Academicos/GitProjects/MarineFoodWebsSmallWorld/Data/FoodWebs")
nets <- readNetwork(dn[20],"~/Academicos/GitProjects/MarineFoodWebsSmallWorld/Data/FoodWebs")
nets <- readNetwork(dn[21],"~/Academicos/GitProjects/MarineFoodWebsSmallWorld/Data/FoodWebs")
V(nets)[[]]
nets <- readNetwork(dn[22],"~/Academicos/GitProjects/MarineFoodWebsSmallWorld/Data/FoodWebs")
nets <- readNetwork(dn[23],"~/Academicos/GitProjects/MarineFoodWebsSmallWorld/Data/FoodWebs")
nets <- readNetwork(dn,"~/Academicos/GitProjects/MarineFoodWebsSmallWorld/Data/FoodWebs")
nets$CaymanIs_FW <- readNetwork(dn[8],"~/Academicos/GitProjects/MarineFoodWebsSmallWorld/Data/FoodWebs",fhead = FALSE,skipColumn = 0)
nets$Jamaica_FW <- readNetwork(dn[15],"~/Academicos/GitProjects/MarineFoodWebsSmallWorld/Data/FoodWebs",fhead = FALSE,skipColumn = 0)
nets$Cuba_FW <- readNetwork(dn[11],"~/Academicos/GitProjects/MarineFoodWebsSmallWorld/Data/FoodWebs",fhead = FALSE,skipColumn = 0)

calcTopologicalIndices(nets)
netData <- nets
devtools::use_data(netData,overwrite = TRUE)

# test reading .dat files
#

dn <- list.files("~/Dropbox/Projects/NetworksAsUnifyingPrinciple/Data",pattern = "^.*\\.dat$")
nets <- readNetwork(dn,"~/Dropbox/Projects/NetworksAsUnifyingPrinciple/Data",fhead=FALSE)
calcTopologicalIndices(list(netData$CaymanIs_FW,nets$cayman_islands))
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
                           3 -+ 6,5 -+8, 8-+ 9, simplify = FALSE)
plotTrophLevel(g,vertexLabel = TRUE,vertexSizeFactor = 20,modules = FALSE)
calcIncoherence(g)
calcModularitySWnessZScore(g,100,0.1,paralell = FALSE)
dg <- components(g)
m<-cluster_spinglass(g)
m$membership

#
# Testing Curve Ball algorithm
#
g <-   graph_from_literal( 1 -+ 4 -+ 7,2 -+ 5 -+ 9, 4+-5,
                           3 -+ 6 -+ 8,5 -+8, simplify = FALSE)
plotTrophLevel(g,vertexLabel = TRUE,vertexSizeFactor = 20,modules = TRUE)
g
gg <- curveBall(g,10)
plotTrophLevel(gg[[1]],vertexLabel = TRUE,vertexSizeFactor = 10,modules = TRUE)
plotTrophLevel(gg[[2]],vertexLabel = TRUE,vertexSizeFactor = 10,modules=TRUE)
plotTrophLevel(gg[[3]],vertexLabel = TRUE,vertexSizeFactor = 10,modules = TRUE)
plotTrophLevel(gg[[10]],vertexLabel = TRUE,vertexSizeFactor = 20,modules = TRUE)
lapply(gg, function(g) components(g)$no)
calcTopologicalIndices(gg)

gg <- curveBall(netData[[1]],10)
lapply(gg,count_components)

# generalization = in degree
# vulnerability  = out degree

dn <- list.files("~/Dropbox/Projects/NetworksAsUnifyingPrinciple/Data",pattern = "^PotterCove.*\\.txt$")
gt <- readMultiplex("~/Dropbox/Projects/NetworksAsUnifyingPrinciple/Data/PotterCove_Multiplex.txt",format="GLV")

E(gt)[E(gt)$type=="Trophic"]

plotTrophLevel(gt,modules=T)



#
#

fileName <- c(system.file("extdata",  package = "EcoNetwork"))
dn <- list.files("inst/extdata",pattern = "^Kefi2015.*\\.txt$")
g <- readNetwork(dn,"inst/extdata", skipColumn = 2)
types <- c("Competitive","Mutualistic","Trophic")
gt <- readMultiplex(dn,types,"inst/extdata", skipColumn = 2)
g[[1]]

gt[]

dg <- degree(g[[1]])
class(dg)
dg2 <- degree(g[[2]])
dg3 <- degree(g[[3]])

sort(degree(gt)) == sort(dg+dg2+dg3)
require(NetIndices)
require(RColorBrewer)
require(igraph)
plotTrophLevel(g[[1]],modules = T)
plotTrophLevel(g[[2]],modules = T)
plotTrophLevel(g[[3]],modules = T)
plotTrophLevel(gt,modules = T)


g <-   graph_from_literal( 2 -+ 1 +-3,4 -+ 1, 4-+4, 3+-3, 5-+5, 4-+6-+2, 2+-5-+3, simplify = FALSE)
E(g)$type <- "Trophic"


pltMat <- plotTrophLevel(g,vertexLabel = TRUE,vertexSizeFactor = 20,modules = TRUE,tk=TRUE)
pltMat[6,2] <- 1
pltMat[6,1] <- 400
pltMat[4,1] <- 20
pltMat[4,2] <- 1
pltMat[5,1] <- 250
pltMat[5,2] <- 70
pltMat[3,1] <- 10
pltMat[3,2] <- 220
pltMat[1,1] <- 450
pltMat[1,2] <- 339
pltMat[2,1] <- 150
pltMat[2,2] <- 376
plotTrophLevel(g,vertexLabel = TRUE,vertexSizeFactor = 20,lMat=pltMat)
