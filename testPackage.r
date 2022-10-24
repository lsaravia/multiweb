# Test
require(igraph)
fileName <- c(system.file("extdata",  package = "multiweb"))
fileName <- system.file("extdata", "BarentsBoreal_FW.csv", package = "multiweb")
g <- readNetwork(fileName)
V(g)$weigth <-  runif(vcount(g))
plotTrophLevel(g,vertexLabel = FALSE,vertexSizeFactor = 3,modules = TRUE,maxTL=4,weights=NULL, frame.color="black")
tk <- plotTrophLevel(g,vertexLabel = FALSE,vertexSizeFactor = 3,modules = TRUE, tk=TRUE)
plotTrophLevel(g,vertexLabel = FALSE,vertexSizeFactor = 3,modules = TRUE, lMat=tk)

dn <- list.files("inst/extdata",pattern = "^.*\\.csv$")
g<- readNetwork(dn,"inst/extdata")
names(g)
E(g[[1]])[which_loop(g[[1]])]

plotTrophLevel(g[[1]])

require(igraph)
require(NetIndices)
g <-   graph_from_literal( 1 -+ 4 -+ 7,2 -+ 5 -+7, 3-+6-+7, 7-+7, 4+-3, simplify = FALSE)
calc_topological_indices(g)
plotTrophLevel(g,vertexLabel = TRUE,vertexSizeFactor = 20,modules = TRUE,maxTL=4,weights=NULL)
plotTrophLevel(g,vertexLabel = TRUE,vertexSizeFactor = 20,modules = TRUE,maxTL=4,weights=NULL, bpal=viridisLite::viridis(11))

# Test tk and use of lMat
tk <- plotTrophLevel(g,vertexLabel = TRUE,vertexSizeFactor = 20,modules = TRUE, tk=TRUE,vertex.label.color="white")
tk [,2] <- seq(0,100, length.out = 7)
plotTrophLevel(g,vertexLabel = TRUE,vertexSizeFactor = 20,modules = TRUE, lMat=tk,maxTL=2)



# add weights
#
E(g)$weight <- c(2,2.5,0.1,0.2,1,5,0.1,2)
plotTrophLevel(g,vertexLabel = TRUE,vertexSizeFactor = 20,modules = TRUE,maxTL=4,weights=NA,edge.width=5)

plotTrophLevel(g,vertexLabel = TRUE,vertexSizeFactor = 20,modules = TRUE,maxTL=4,weights=NULL,edge.width=5)


# Test incoherence
#
calcIncoherence(g)

g <-   graph_from_literal( 1 -+ 4 -+ 7,2 -+ 5 -+7, 3-+6-+7, 4+-5, simplify = FALSE)
calcTopologicalIndices(g)
plotTrophLevel(g,vertexLabel = TRUE,vertexSizeFactor = 20, modules = TRUE)
nullg <- generateERbasal(g,10)
calcIncoherence(nullg,ncores=0)

calcModularitySWnessZScore(g,nullg,0.1)
calc_swness_zscore(g,nullg,0.1,ncores = 0)

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
calcIncoherence(nets)
#netData <- nets
#devtools::use_data(netData,overwrite = TRUE)


# test reading .dat files from Johnson 2017
#

dn <- list.files("~/Dropbox/Projects/NetworksAsUnifyingPrinciple/Data",pattern = "^.*\\.dat$")
nets <- readNetwork(dn,"~/Dropbox/Projects/NetworksAsUnifyingPrinciple/Data",fhead=FALSE)
calcTopologicalIndices(list(netData$CaymanIs_FW,nets$cayman_islands))
require(RColorBrewer)
plotTrophLevel(nets$cayman_islands)
plotTrophLevel(netData$CaymanIs_FW)


#
nets <- readNetwork("Meta.dat","~/Dropbox/Projects/StabilityConstraints/Data",edgeListFormat = 2)



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
calcModularitySWnessZScore(g,generateERbasal(g,10),0.1)


# Completely coherent Q=0 network
#
g <-   graph_from_literal( 1 -+ 4 -+ 7,2 -+ 5 , 2-+4, 5-+7, 3-+6 -+8,
                           5 -+8, simplify = FALSE)
plotTrophLevel(g,vertexLabel = TRUE,vertexSizeFactor = 20,modules = TRUE)
calcTopologicalIndices(g)
calcIncoherence(g)
calcModularitySWnessZScore(g,generateERbasal(g,100),0.1,ncores = 0)
dg <- components(g)
m<-cluster_spinglass(g)
m$membership

#
# Testing Curve Ball algorithm
#
g <-   graph_from_literal( 2 -+ 1 +-3,4 -+ 1, 3+-3, 4-+6-+2, 2+-5-+3, simplify = FALSE)
calcTopologicalIndices(g)
plotTrophLevel(g,vertexLabel = TRUE,vertexSizeFactor = 20,modules = TRUE)
gg <- curveBall(g,10)
gg <- curve_ball(g,1000)
calcTopologicalIndices(gg,ncores=48)

plotTrophLevel(gg[[1]],vertexLabel = TRUE,vertexSizeFactor = 20,modules = TRUE)
plotTrophLevel(gg[[2]],vertexLabel = TRUE,vertexSizeFactor = 10,modules=TRUE)
plotTrophLevel(gg[[3]],vertexLabel = TRUE,vertexSizeFactor = 10,modules = TRUE)
plotTrophLevel(gg[[10]],vertexLabel = TRUE,vertexSizeFactor = 20,modules = TRUE)


#
# Add weight
#
E(g)$weight <- sample(c(.1,.2,.8,.9),gsize(g),replace=TRUE)
m <- get.adjacency(g,sparse=FALSE,attr="weight")
gg <- curveBall(g,2,istrength =TRUE)
m1 <- get.adjacency(gg[[1]],sparse=FALSE)
m
m1

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

dn <- list.files("~/Academicos/InconclusedProjects/NetworksAsUnifyingPrinciple/Data",pattern = "^PotterCove.*\\.txt$")
gt <- readMultiplex("~/Academicos/InconclusedProjects/NetworksAsUnifyingPrinciple/Data/PotterCove_Multiplex.txt",format="GLV")

# Number of interaction by type
sapply(gt,ecount)
sapply(gt,vcount)

# Proportion of interactions by type
sapply(gt,ecount)/sum(sapply(gt,ecount))

# Connectivity
sapply(gt,ecount)/(sapply(gt,vcount)*sapply(gt,vcount))

# Multiplex QSS
calc_QSS(gt)
#   QSS    MEing
#    0    25.87422
calc_QSS(gt[[3]])
#   QSS    MEing
#    0    2.099394

glv <- toGLVadjMat(gt)

sum(glv==-1)   # Competitive + trophic


# Meta-web assembly models package
#
require(meweasmo)


pin  <- calcPropInteractionsGLVadjMat(glv, rep(1,times=nrow(glv)))

# Read the file directly
mtm <- as.matrix(mtm[,2:92])
#
mtm <- read.delim("~/Dropbox/Projects/NetworksAsUnifyingPrinciple/Data/PotterCove_Multiplex.txt")
pin1  <- calcPropInteractionsGLVadjMat(mtm, rep(1,times=nrow(glv)))

# The proportion of types of interactions must be equal
#
pin == pin1


plotTrophLevel(gt[[3]],modules = T)


# Read multiple interaction network in different layers
#

fileName <- c(system.file("extdata",  package = "multiweb"))
dn <- list.files(fileName,pattern = "^Kefi2015.*\\.txt$")
g <- readNetwork(dn,fileName, skipColumn = 2)
gt <- fromIgraphToMgraph(g,c("Negative","Positive","Antagonistic"))

types <- c("Competitive","Mutualistic","Trophic")
gt <- readMultiplex(dn,types,fileName, skipColumn = 2)
tp <- calc_topological_indices(gt)
# how many interactions from each type
#
sapply(gt,ecount)/sum(sapply(gt,ecount))

# Number of species
sapply(gt,vcount)

# Total Connectivity
sum(sapply(gt,ecount))/(106*106)


# Calc multi type interaction QSS
#
calc_QSS(gt)
# QSS    MEing
# 1   0 26.04826

calc_QSS(g[[3]])               # QSS with trophic interactions only should be more stable
#  QSS    MEing
# 1   0 3.434581



#
# Select the first type of interaction
#
nm <- names(gt)
gt[[nm[1]]]

glv <- toGLVadjMat(gt)

sum(glv==-1)   # Competitive + trophic

#
# Simple test
#
g <-   graph_from_literal( 2 -+ 1 +-3,4 -+ 1, 4-+4, 3+-3, 5-+5, 4-+6-+2, 2+-5-+3, simplify = FALSE)
E(g)$type <- "Trophic"
c <- cluster_infomap(as.undirected(g))
plot_troph_level(g,vertexLabel = TRUE,vertexSizeFactor = 20,vertexSizeMin = 12, modules = TRUE, community_obj = c)
plot_troph_level(g,tk=FALSE,vertexLabel = TRUE,vertexSizeFactor = 20,vertexSizeMin = 12, modules = TRUE,vertex.label.color="white",edge.width = 0.3)
plot_troph_level(g,modules =FALSE ,tk=TRUE, vertexLabel = TRUE,vertexSizeFactor = 20)

g1 <-   graph_from_literal(1,2,3, 4 -+ 5, 5-+4,6, simplify = FALSE)
E(g1)$type <- "Competitive"

g2 <-   graph_from_literal(1 -+ 2,2,3,4,5, 6 -+ 3, 3-+6, simplify = FALSE)
E(g2)$type <- "Mutualistic"
g2
#
# Builds an mgraph object
#
mg <- fromIgraphToMgraph(list(g1,g2,g),c("Competitive", "Mutualistic", "Trophic"))

#
# Automatically uses the function toGLVadjMat to calculate the Generalized Lotka-Volterra interaction matrix
# and then calculates the QSS
#
set.seed(423)
calc_QSS(list(g,g1,g2))
set.seed(423)
calc_QSS(list(g,g1,g2),returnRaw = TRUE)


calc_QSS(mg)

#
# Add weights
#
set.seed(1231)
E(g)$weight <- sample(c(.1,.2,.3,.9),gsize(g),replace=TRUE)
E(g1)$weight <- sample(c(.1,.2,.8,.9),gsize(g1),replace=TRUE)
E(g2)$weight <- sample(c(.1,.2,.8,.9),gsize(g2),replace=TRUE)
mg <- fromIgraphToMgraph(list(g1,g2,g),c("Competitive", "Mutualistic", "Trophic"))
toGLVadjMat(mg,istrength=TRUE)
calc_QSS(list(g,g1,g2))calc_QSS(list(g,g1,g2),istrength = TRUE)
calc_QSS(mg)
calc_QSS(mg,istrength = TRUE)
calc_QSS(mg,istrength = TRUE,nsim=1)
calc_QSS(mg,istrength = FALSE,nsim=1)

# Calc QSS with mean & maximum weigth
#
set.seed(3231)
calc_QSS(g,istrength = FALSE)
calc_QSS(g,istrength = TRUE)
calc_QSS(g,istrength = TRUE,nsim = 1)
calc_QSS(g,istrength = FALSE,nsim = 1)

g1 <-  g
E(g1)$weight <- max(E(g)$weight)
calc_QSS(g1,istrength = TRUE)
calc_QSS(g1,istrength = TRUE,nsim=1)

E(g1)$weight <- mean(E(g)$weight)
calc_QSS(g1,istrength = TRUE)

#
require(meweasmo)

calcPropInteractionsGLVadjMat(glv, rep(1,times=nrow(glv)))



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

#
# Test calc_QSS
#
set.seed(342)
calc_QSS(g,istrength = TRUE, nsim=10)
set.seed(342)
calc_QSS(g,istrength = TRUE, nsim=10,returnRaw = TRUE)

glvI <- fromIgraphToMgraph(list(make_empty_graph(n=vcount(g)),make_empty_graph(n=vcount(g)),g) ,c("empty","empty","Antagonistic"))
glv <- toGLVadjMat(glvI,c("empty","empty","Antagonistic"),istrength = TRUE)   #
mean(glv[glv>0])
glv[glv>0] <- glv[glv>0]*0.1
plotTrophLevel(g,vertexLabel = TRUE,vertexSizeFactor = 5,modules = TRUE)

null <- curveBall(g,100)
calc_weighted_topological_indices(null[[1]])
calc_topological_indices(null[[2]])

null <- curveBall(g,100,istrength = TRUE)
calc_weighted_topological_indices(null[[1]])
calc_weighted_topological_indices(null[[2]])

calc_modularity_swness_zscore(g,null)

# Plot a network with two different components as modules
#
g <-   graph_from_literal( 2 -+ 1 +-3,4 -+ 1, 3+-3, 5-+5, 4-+6-+2, 2+-5-+3, simplify = FALSE)
g1 <-   graph_from_literal( 7 -+ 8 +-9,10 -+ 8, 10-+10, 9+-9, 11-+11, 10-+13-+12, 12+-11-+9, simplify = FALSE)

g2 <- g+g1
plotTrophLevel(g2,vertexLabel = TRUE,vertexSizeFactor = 7,modules = TRUE)

#
# Test calc_weighted_topological_indices
#

E(g)$weight <- sample(c(0.001,.01,.1,.5,.9),gsize(g),replace=TRUE)
calc_weighted_topological_indices(g)
E(g)$weight <- rep(.01,times=gsize(g))
calc_weighted_topological_indices(g)
calc_topological_indices(g)


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
m <- matrix(0,nrow=5,ncol=5)
m[1,2] <- m[1,3] <- m[3,4]<- .2
m[2,1] <- m[3,1] <- m[4,3] <- -2
m[5,4] <- m[4,5] <- 0.1          # Mutualistic
m[1,1] <- -0.01                    # Cannibalistic
g <- fromGLVadjToIgraph(m,c(1,1,1,1,0))
g <- fromGLVadjToIgraph(m,c(0,1,1,1,1))

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

sapply(netData,igraph::vcount)

g <- netData[[1]]

tp <- calc_topological_roles(g,nsim=10,ncores=4)

classify_topological_roles(tp,g,plt=TRUE)



#
# Test QSS
#

sapply(netData,vcount)

set.seed(342)
g <- netData[[2]]
g1 <- netData[[18]]
tp <- calc_QSS(list(g,g1),nsim=100,ncores=8)


# with default values QSS is almost always 0
tp
#  QSS    MEing
#1   0 4.591426
#2   0 2.052344

set.seed(342)
tp <- calc_QSS(list(g,g1),nsim=100,ncores=8,negative=-1,positive=1)
tp
#      QSS      MEing
# 1 0.0005 0.24024671
# 2 0.3224 0.07553226

# Test for differences
#
prop.test(x=tp$QSS*10000, n=c(10000,10000))


#
# test calc_modularity
#

nullg <- generateERbasal(netData[[1]],10)
calc_modularity(nullg)

# with weigths
g <- netData[[1]]
names(netData)[1]
calc_modularity(g)

V(g)$weight <-  runif(vcount(g))
calc_modularity(g,weights = NULL)

#
# test QSS with 1 simulation vs mean of QSS with 1000 simulations
#
set.seed(5431)
g <- netData[[1]]
E(g)$weight <- runif(ecount(g))

calc_QSS(g,istrength = FALSE)                  # 1   0 3.106538
calc_QSS(g,istrength = TRUE,ncores =4 )        # 1   0 1.575575
calc_QSS(g,istrength = FALSE,nsim=1)           # 1   0 2.461278
calc_QSS(g,istrength = TRUE,nsim=1)            # 1   0 1.421724

raw <-  calc_QSS(g,istrength = TRUE,ncores =4 ,returnRaw = TRUE)        # 1   0 1.637529
require(ggplot2)
ggplot(raw, aes(maxre)) + geom_density() + theme_bw() + geom_vline(xintercept = calc_QSS(g,istrength = TRUE,nsim=1)$MEing, linetype="dashed")

#
# using lognormal interaction strengths
#

E(g)$weight <- rlnorm(ecount(g))
calc_QSS(g,istrength = FALSE)                  # 1   0 3.093221
calc_QSS(g,istrength = TRUE,ncores =4 )        # 1   0 4.577607
calc_QSS(g,istrength = FALSE,nsim=1)           # 1   0 2.461278
calc_QSS(g,istrength = TRUE,nsim=1)            # 1   0 4.01689
raw <-  calc_QSS(g,istrength = TRUE,ncores =4 ,returnRaw = TRUE)
ggplot(raw, aes(maxre)) + geom_density() + theme_bw() + geom_vline(xintercept = calc_QSS(g,istrength = TRUE,nsim=1)$MEing, linetype="dashed")


#
# Test calc_QSS_extinction_dif
#

g <- netData[[1]]
E(g)$weight <-  runif(vcount(g))
E(g)$w <-  runif(vcount(g))
calc_QSS(g,istrength = FALSE)
calc_QSS(g,istrength = TRUE)

edge_attr_names(g)

calc_QSS_extinction_dif(g,V(g)$name[1:3],nsim=10,istrength = FALSE)

calc_QSS_extinction_dif(g,V(g)$name[1:3],nsim=10,istrength = TRUE)

#
# Test calc_QSS_extinction_dif
#
eseq <- calc_QSS_extinctions_seq(g,V(g)$name[1:5],nsim=10,istrength = FALSE)
eseq

eseq <- calc_QSS_extinctions_seq(g,V(g)$name,nsim=100,istrength = TRUE)
eseq



# Test calc_interaction_intensity
#
# get the data frame of interactions
#
g <- netData[[1]]
require(dplyr)
require(igraph)
set.seed(7815)
da <- as_long_data_frame(g) %>% dplyr::select(from:to) %>% mutate(con_mm=rlnorm(n(),5,2),res_mm=con_mm - 30 ,int_dim=sample(c("2D","3D"),n(),replace=TRUE), res_den = -999)
calc_interaction_intensity(da,res_mm,res_den,con_mm,int_dim)

# Using values for res_den
#
da <- as_long_data_frame(g) %>% dplyr::select(from:to) %>%
  mutate(con_mm=rlnorm(n(),5,2),res_mm=con_mm - 30 ,int_dim=sample(c("2D","3D"),n(),replace=TRUE), res_den = runif(n(),1e-5,1))

calc_interaction_intensity(da,res_mm,res_den,con_mm,int_dim)


