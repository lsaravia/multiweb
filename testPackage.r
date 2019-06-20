# Test
require(igraph)
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
print(calc_topological_indices(g))
plotTrophLevel(g,vertexLabel = TRUE,vertexSizeFactor = 20)
calcIncoherence(g)

g <-   graph_from_literal( 1 -+ 4 -+ 7,2 -+ 5 -+7, 3-+6-+7, 4+-5, simplify = FALSE)
calcTopologicalIndices(g)
plotTrophLevel(g,vertexLabel = TRUE,vertexSizeFactor = 20, modules = TRUE)
calcIncoherence(g)
nullg <- generateERbasal(g,10)
calcModularitySWnessZScore(g,nullg,0.1,ncores = NULL)


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

# test reading .dat files from Johnson 2017
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
calcTopologicalIndices(nets)
calcIncoherence(nets)
mean(degree(nets$benguela,mode="out"))
calcIncoherence(nets$bridge)

mean(degree(nets$bridge))
degree()

# Test modularity & incoherence Q
# 1 =Omnivory link 5+-4
#
g <-   graph_from_literal( 1 -+ 4 -+ 7,2 -+ 5 -+ 9, 4+-5,
                           3 -+ 6 -+ 8,5 -+8, simplify = FALSE)
plotTrophLevel(g,vertexLabel = TRUE,vertexSizeFactor = 20,modules = TRUE)
title(main="9 species")
calcTopologicalIndices(g)
calcIncoherence(g)
calcModularitySWnessZScore(g,100,0.1)


# Completely coherent Q=0 network
#
g <-   graph_from_literal( 1 -+ 4 -+ 7,2 -+ 5 , 2-+4, 5-+7, 3-+6 -+8,
                           5 -+8, simplify = FALSE)
plotTrophLevel(g,vertexLabel = TRUE,vertexSizeFactor = 20,modules = TRUE)
calcTopologicalIndices(g)
calcIncoherence(g)
calcModularitySWnessZScore(g,100,0.1,ncores = 0)
dg <- components(g)
m<-cluster_spinglass(g)
m$membership

#
# Testing Curve Ball algorithm
#
g <-   graph_from_literal( 1 -+ 4 -+ 7,2 -+ 5 -+ 9, 4+-5,
                           3 -+ 6 -+ 8,5 -+8, simplify = FALSE)
plotTrophLevel(g,vertexLabel = TRUE,vertexSizeFactor = 20,modules = TRUE)
gg <- curveBall(g,100)
plotTrophLevel(gg[[1]],vertexLabel = TRUE,vertexSizeFactor = 20,modules = TRUE)
plotTrophLevel(gg[[2]],vertexLabel = TRUE,vertexSizeFactor = 10,modules=TRUE)
plotTrophLevel(gg[[3]],vertexLabel = TRUE,vertexSizeFactor = 10,modules = TRUE)
plotTrophLevel(gg[[10]],vertexLabel = TRUE,vertexSizeFactor = 20,modules = TRUE)

# Test SW zscores
#
any(sapply(gg, function(g) components(g)$no>1))
calcModularitySWnessZScore(g,gg,0.01)

# Test Parallel
#
p <- calcTopologicalIndices(gg,4)
p <- calcIncoherence(gg,4)

# Calculate 1000 randomizations of a real network
#
gg <- curveBall(netData[[1]],100)
sapply(gg, function(g) components(g)$no)

lapply(gg,count_components)!=1
tp <- calcTopologicalIndices(gg)
tp1<- calcTopologicalIndices(gg,4)
any(tp!=tp1)
tp <- calcIncoherence(gg)
tp1 <- calcIncoherence(gg,ncores=4)
any(tp!=tp1)

#
# Duplicate the same network to test the function
#
gg <- rep(list(netData[[1]]),100)
p <- calcIncoherence(gg)
p1 <- calcIncoherence(gg,ncores=4)
any(p!=p1)


require(ggplot2)
ggplot(tp, aes(Omnivory)) + geom_density() + theme_bw()

require(ggplot2)
ggplot(tp, aes(Q)) + geom_density() + theme_bw()
ggplot(tp, aes(mTI)) + geom_density() + theme_bw()

tp1<- calcIncoherence(netData[[1]])
tp1$y <- 0
ggplot(tp, aes(Q)) + geom_density() + theme_bw() + geom_point(data=tp1,aes(Q,y))


# generalization = in degree
# vulnerability  = out degree

dn <- list.files("~/Dropbox/Projects/NetworksAsUnifyingPrinciple/Data",pattern = "^PotterCove.*\\.txt$")
gt <- readMultiplex("~/Dropbox/Projects/NetworksAsUnifyingPrinciple/Data/PotterCove_Multiplex.txt",format="GLV")

# Number of interaction by type
sapply(gt,ecount)
sapply(gt,vcount)

# Proportion of interactions by type
sapply(gt,ecount)/sum(sapply(gt,ecount))

# Connectivity
sapply(gt,ecount)/(sapply(gt,vcount)*sapply(gt,vcount))


glv <- toGLVadjMat(gt)

sum(glv==-1)   # Competitive + trophic

require(MetaWebAssemblyModels)


pin  <- calcPropInteractionsGLVadjMat(glv, rep(1,times=nrow(glv)))

# Read the file directly
#
mtm <- read.delim("~/Dropbox/Projects/NetworksAsUnifyingPrinciple/Data/PotterCove_Multiplex.txt")
mtm <- as.matrix(mtm[,2:92])
pin1  <- calcPropInteractionsGLVadjMat(mtm, rep(1,times=nrow(glv)))

# The proportion of types of interactions must be equal
#
pin == pin1


plotTrophLevel(gt[[3]],modules = T)


# Read multiple interaction network in different layers
#

fileName <- c(system.file("extdata",  package = "EcoNetwork"))
dn <- list.files("inst/extdata",pattern = "^Kefi2015.*\\.txt$")
g <- readNetwork(dn,"inst/extdata", skipColumn = 2)
gt <- igraph2mgraph(g,c("Negative","Positive","Antagonistic"))

types <- c("Competitive","Mutualistic","Trophic")
gt <- readMultiplex(dn,types,"inst/extdata", skipColumn = 2)

# how many interactions from each type
#
sapply(gt,ecount)/sum(sapply(gt,ecount))

# Number of species
sapply(gt,vcount)

# Total Connectivity
sum(sapply(gt,ecount))/(106*106)

#
# Select the first type of interaction
#
nm <- names(gt)
gt[[nm[1]]]

glv <- toGLVadjMat(gt)

sum(glv==-1)   # Competitive + trophic

require(MetaWebAssemblyModels)

calcPropInteractionsGLVadjMat(glv, rep(1,times=nrow(glv)))

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

#
# Test network with weights
#
fileName <- "../../Collaborations/NetworkGolfoSanJorge/Data/TTSinPescaInteraccion.dat"
g <- readNetwork(fileName,edgeListFormat=2 )
V(g)[nei(V(g)[name=="POM"],"out")]
V(g)[nei(V(g)[name=="POM"],"in")]
g
g[]
# Get the weights of a predator
g[,"Acanthistius patachonicus"]
# Get the adjacency matrix with weights
as_adjacency_matrix(g,attr="weight",sparse = FALSE)

edge_attr(g,"weight")

glvI <- fromIgraphToMgraph(list(make_empty_graph(n=vcount(g)),make_empty_graph(n=vcount(g)),g) ,c("empty","empty","Antagonistic"))
glv <- toGLVadjMat(glvI,c("empty","empty","Antagonistic"),istrength = TRUE)   #
mean(glv[glv>0])
glv[glv>0] <- glv[glv>0]*0.1
plotTrophLevel(g,vertexLabel = TRUE,vertexSizeFactor = 5,modules = TRUE)

null <- curveBall(g,100)
calc_modularity_swness_zscore(g,null)

# Plot a network with two different components as modules
#
g <-   graph_from_literal( 2 -+ 1 +-3,4 -+ 1, 3+-3, 5-+5, 4-+6-+2, 2+-5-+3, simplify = FALSE)
g1 <-   graph_from_literal( 7 -+ 8 +-9,10 -+ 8, 10-+10, 9+-9, 11-+11, 10-+13-+12, 12+-11-+9, simplify = FALSE)

g2 <- g+g1
plotTrophLevel(g2,vertexLabel = TRUE,vertexSizeFactor = 7,modules = TRUE)

# Test quantitative connectance
#
#
# Equal fluxes
m <- matrix(0,nrow=4,ncol=4)
m[1,2] <- m[1,3] <- m[3,4]<- 2
m[2,1] <- m[3,1] <- m[4,3] <- -2
calc_quantitative_connectance(m, c(1,1,1,1))
calcQuantitativeConnectance(m, c(1,1,1,1))
g <- fromGLVadjToIgraph(m,c(1,1,1,1))
print(calcTopologicalIndices(g))
print(calc_topological_indices(g))
plotTrophLevel(g,vertexLabel = T)
# Unequal fluxes
#
m <- matrix(0,nrow=4,ncol=4)
m[1,2] <- m[1,3] <- m[3,4]<- .2
m[2,1] <- m[3,1] <- m[4,3] <- -2
calc_quantitative_connectance(m, c(1,1,1,1))

# Examples from  Ulanowicz, R.E. & Wolff, W.F. (1991). Ecosystem flow networks: Loaded dice? Math. Biosci., 103, 45–68
#
m <- matrix(0,nrow=4,ncol=4)
m[1,2] <- m[2,3] <- m[3,4] <- m[4,1] <- 75
m[1,3] <- m[3,1] <- m[2,4] <- m[4,2] <- 75
calc_quantitative_connectance(m, c(1,1,1,1)) # 0.5, 2
sum(m>0)/(nrow(m)^2)


m <- matrix(0,nrow=4,ncol=4)
m[1,2] <- m[2,3] <- m[3,4] <- m[4,1] <- 50
m[1,3] <- m[3,1] <- m[2,4] <- m[4,2] <- 100
print(calc_quantitative_connectance(m, c(1,1,1,1))) # =0.4724704, 1.889882
sum(m>0)/(nrow(m)^2)

m <- matrix(0,nrow=4,ncol=4)
m[1,2] <-  560
m[1,3] <- m[2,3] <- m[3,4] <- m[4,1] <- m[3,1] <- m[2,4] <- m[4,2] <- 5.66
print(calc_quantitative_connectance(m, c(1,1,1,1))) # =0.4411018, 1.096489
sum(m>0)/(nrow(m)^2)

#
# Test Topological roles
#

sapply(netData,vcount)

g <- netData[[2]]

tp <- calc_topological_roles(g,nsim=100,ncores=4)

classify_topological_roles(tp,g,plt=TRUE)



#
# Test QSS
#

sapply(netData,vcount)

set.seed(342)
g <- netData[[2]]
g1 <- netData[[18]]
tp <- calc_QSS(list(g,g1),nsim=10000,ncores=48)

prop.test(x=tp$QSS*10000, n=c(10000,10000))

tp
#      QSS      MEing
# 1 0.0005 0.24024671
# 2 0.3224 0.07553226


#
# test calc_modularity
#

nullg <- generateERbasal(netData[[1]],10)
calc_modularity(nullg)


#
#
#
nullg <- generateERbasal(netData[[1]],10)
calc_modularity(nullg)
