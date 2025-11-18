# Test
require(igraph)
fileName <- c(system.file("extdata",  package = "multiweb"))
fileName <- system.file("extdata", "BarentsBoreal_FW.csv", package = "multiweb")
g <- readNetwork(fileName)
V(g)$weigth <-  runif(vcount(g))
plotTrophLevel(g,vertexLabel = FALSE,vertexSizeFactor = 3,modules = TRUE,maxTL=4,weights=NULL, frame.color="black")
tk <- plotTrophLevel(g,vertexLabel = FALSE,vertexSizeFactor = 3,modules = TRUE, tk=TRUE)
plotTrophLevel(g,vertexLabel = FALSE,vertexSizeFactor = 3,modules = TRUE, lMat=tk,tk=TRUE)

#
dn <- list.files("inst/extdata",pattern = "^.*\\.csv$")
g<- readNetwork(dn,"inst/extdata")
names(g)

plotTrophLevel(g[[1]])

require(igraph)
require(NetIndices)
g <-   graph_from_literal( 1 -+ 4 -+ 7,2 -+ 5 -+7, 3-+6-+7, 7-+7, 4+-3, 2-+7, simplify = FALSE)
calc_topological_indices(g)
plotTrophLevel(g,vertexLabel = TRUE,vertexSizeFactor = 20,modules = TRUE,maxTL=4,weights=NULL)
plotTrophLevel(g,vertexLabel = TRUE,vertexSizeFactor = 20,modules = TRUE,maxTL=4,weights=NULL, bpal=viridisLite::viridis(11))

# Test tk and use of lMat
tk <- plotTrophLevel(g,vertexLabel = TRUE,vertexSizeFactor = 20,modules = TRUE, tk=TRUE,vertex.label.color="white")
tk [1,2] <- 100 # seq(0,100, length.out = 7)
plotTrophLevel(g,vertexLabel = TRUE,vertexSizeFactor = 20,modules = TRUE, lMat=tk,maxTL=3)

V(g)$weight <- c(100,2.5,0.1,0.2,1,5,0.1)
plot_troph_level_ggraph(g,node_weights=NA)
plot_troph_level_ggplot(g,node_weights=NA)

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
plotTrophLevel(nullg[[10]],vertexLabel = TRUE,vertexSizeFactor = 20, modules = TRUE)
plot_troph_level_ggraph(nullg[[1]])
calcIncoherence(g,ncores=0)
calcIncoherence(nullg,ncores=0)

nullg <- generate_niche(7,0.1428571,10)

calcModularitySWnessZScore(g,nullg,0.1)
calc_swness_zscore(g,nullg,0.1,ncores = 0)

# test reading .csv files using a folder with multiple csv files
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


# test reading .dat files from Johnson 2017
#
dn <- list.files("~/Dropbox/Projects/NetworksAsUnifyingPrinciple/Data",pattern = "^.*\\.dat$")
nets <- readNetwork(dn,"~/Dropbox/Projects/NetworksAsUnifyingPrinciple/Data",fhead=FALSE)
calcTopologicalIndices(list(netData$CaymanIs_FW,nets$cayman_islands))


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
gg <- curveBall(g,10,istrength =TRUE)
m1 <- get.adjacency(gg[[1]],sparse=FALSE,attr="weight")
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



# Test reading multiple interaction networks
#
fileName <- c(system.file("extdata",  package = "multiweb"))
dn <- list.files(fileName,pattern = "^Kefi2015.*\\.txt$")
gt <- readMultiplex(dn,types,"inst/extdata", skipColumn = 2)

# Number of interaction by type
sapply(gt,ecount)
sapply(gt,vcount)

# Proportion of interactions by type
sapply(gt,ecount)/sum(sapply(gt,ecount))

# Connectivity
sapply(gt,ecount)/(sapply(gt,vcount)*sapply(gt,vcount))

calc_centrality(gt)
# Multiplex QSS
calc_QSS(gt)
#   QSS    MEing
#    0    25.87422
calc_QSS(gt[[3]],selfDamping = -10)
#   QSS            MEing
# 0.7 - 0.6     -0.04,-0.03


# Meta-web assembly models package
#
require(meweasmo)

pin  <- calcPropInteractionsGLVadjMat(glv, rep(1,times=nrow(glv)))

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
calc_QSS(gt, selfDamping = 0)
# QSS    MEing
# 1   0 26.04826

calc_QSS(g[[3]],selfDamping = -10,nsim=10)               # QSS with trophic interactions only should be more stable
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
plotTrophLevel(g1,vertexLabel = TRUE,vertexSizeFactor = 7,modules = TRUE)
plot_troph_level_ggplot(g,modules = TRUE, arrow_size = 0.08)
plot_troph_level_ggplot(g1,modules = TRUE)
g <- netData[[1]]

g2 <- g+g1
plotTrophLevel(g2,vertexLabel = TRUE,vertexSizeFactor = 7,modules = TRUE)
plotTrophLevel(g2,vertexLabel = TRUE,vertexSizeFactor = 10,modules = TRUE)
plot_troph_level_ggplot(g2,modules = TRUE)

#
# Test calc_weighted_topological_indices
#

E(g)$weight <- sample(c(0.001,.01,.1,.5,.9),gsize(g),replace=TRUE)
plot_troph_level_ggplot(g,modules = TRUE)
run_infomap(g, return_df = TRUE)
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
g <- netData[[23]]

tp <- calc_topological_roles(g,nsim=10,ncores=4)
classify_topological_roles(tp,g,plt=TRUE)

m <- run_infomap(g)
length(m)
m <- igraph::cluster_walktrap(g)
length(m)
igraph::membership(m)
m <- igraph::cluster_spinglass(g)

tp1 <- multiweb::calc_topological_roles(g,nsim=10,ncores=4,community=m)

multiweb::classify_topological_roles(tp1,g,community=m)



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
E(g)$weight <-  runif(ecount(g))
E(g)$w <-  runif(ecount(g))
calc_QSS(g,istrength = FALSE)
calc_QSS(g,istrength = TRUE)

edge_attr_names(g)

tic()
dif <- calc_QSS_extinction_dif(g,V(g)$name[1:3],nsim=20,istrength = FALSE)
toc()

ggplot(dif, aes(difQSS,fill=Deleted)) + geom_density( alpha=0.3 ) + theme_bw() + geom_vline(xintercept = 0, linetype="dashed") +
  facet_wrap(~Deleted, ncol = 1)

dif1 <- calc_QSS_extinction_dif(g,V(g)$name[1:3],nsim=10,istrength = TRUE)


calc_QSS_extinction_dif_grp(g,V(g)$name[1:3],nsim=10,istrength = TRUE)

#
# Test calc_QSS_extinction_dif
#
eseq <- calc_QSS_extinctions_seq(g,V(g)$name[1:5],nsim=10,istrength = FALSE)
eseq

eseq <- calc_QSS_extinctions_seq(g,V(g)$name,nsim=100,istrength = TRUE)
eseq

g_del <- calc_QSS_extinction_dif_grp(g, bygroup,nsim,ncores,istrength)


# Test calc_interaction_intensity
#
# get the data frame of interactions
#
g <- netData[[1]]
require(dplyr)
require(igraph)
set.seed(7815)
da <- as_long_data_frame(g) %>% dplyr::select(from:to) %>% mutate(con_mm=rlnorm(n(),5,2),res_mm=con_mm - 30 ,int_dim=sample(c("2D","3D"),n(),replace=TRUE), res_den = -999)
calc_interaction_intensity(da,res_mm,res_den,con_mm,int_dim,10)
#
# Test the nsims parameter
#
i1 <- calc_interaction_intensity(da[1,],res_mm,res_den,con_mm,int_dim,1000)
i <- calc_interaction_intensity(da[1,],res_mm,res_den,con_mm,int_dim,1)

require(ggplot2)
ggplot(i1,aes(qRC)) + geom_histogram()
ggplot(i1,aes(qRC)) + geom_density() + scale_x_log10() + geom_vline(xintercept = i$qRC)


# Using values for res_den
#
da <- as_long_data_frame(g) %>% dplyr::select(from:to) %>%
  mutate(con_mm=rlnorm(n(),5,2),res_mm=con_mm - 30 ,int_dim=sample(c("2D","3D"),n(),replace=TRUE), res_den = runif(n(),1e-5,1))

calc_interaction_intensity(da,res_mm,res_den,con_mm,int_dim)


# Test calc_interaction_intensity2
#
# Define an edge list (consumer, resource, bodymass)
edge_list <- tibble(
  predator = c("A", "A", "B", "C"),
  prey = c("B", "C", "D", "D"),
  predator_mass = c(100, 10, 1, .1)  # Body mass for consumers
)

# Compute interaction strength as an edge list (predator effects only)
interaction_strength_edgelist <- calc_interaction_intensity2(
  edge_list, consumer_n = predator, resource_n = prey, bodymass_n = predator_mass, output_format = "edgelist"
)
print(interaction_strength_edgelist)

# Compute interaction strength as an adjacency matrix
interaction_strength_matrix <- calc_interaction_intensity2(
  edge_list, predator,  prey, predator_mass, output_format = "matrix"
)
print(interaction_strength_matrix)

#
# Niche Model
#

generate_niche(20, 0.1)        # Single adjacency matrix
calc_topological_indices(generate_niche(20, 0.1, nsim=5))


#
# Calculate centrality SVD and Shuffling - Niche model
#
calc_svd_entropy_importance(generate_niche(20, 0.1))
calc_svd_entropy_importance(netData[[29]])
calc_centrality(netData[[29]])
calc_centrality(netData[[29]], centrality_func = page_rank)

A <- generate_shuffled_seq_tol(netData[[29]], weighted = FALSE, shuffle_func = shuffle_network_ws)
A$Metrics
A <- generate_shuffled_seq_tol(netData[[29]], weighted = TRUE)

calc_svd_entropy(A$New_A)
calc_svd_entropy(netData[[29]])

calc_svd_entropy(netData[[19]])
B <- shuffle_network_ws(netData[[19]], weighted = FALSE,delta=10)
gB <- graph_from_adjacency_matrix(B, mode = "directed",)
plot_troph_level_ggplot(netData[[19]])
plot_troph_level_ggplot(gB)
calc_svd_entropy(B)

generate_shuffled_seq(netData[[19]], shuffle_func = shuffle_network_ws, weighted = FALSE)
generate_shuffled_seq_tol(netData[[19]], weighted = FALSE, shuffle_func = shuffle_network_ws)

g <- netData[[23]]
E(g)$weight <-  1
generate_shuffled_seq_tol(g, weighted = TRUE, shuffle_func = shuffle_network_ws)
generate_shuffled_seq(g, shuffle_func = shuffle_network_ws, weighted = TRUE)

#
# InfoMap monolayer
#

g <- netData[[23]]
py_infomap0 <- run_infomap(g, output_dir = ".",return_df = TRUE)
py_infomap0$codelength
py_infomap0 <-py_infomap0$communities
py_infomap0 %>% filter(module == 1)
(py_infomap <- run_infomap(g))
membership(py_infomap)
py_infomap$codelength
ig_infomap <- cluster_infomap(g)
ig_infomap$codelength

E(g)$weight <-  runif(ecount(g),0.1,2)
run_infomap(g)

g <- netData[[23]]
gl <- curve_ball(g)
modl <- calc_modularity(gl,cluster_function = run_infomap)


ggplot(modl,aes(x=Modularity)) + geom_density()
ggplot(modl, aes(x = Modularity)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = py_infomap$modularity, color = "red", linetype = "dashed") +
  labs(title = "Distribution of Modularity under Curveball Randomization",
       x = "Modularity", y = "Density") + theme_bw()
#ggsave("~/Downloads/modularity_density_plot.png", width=8,height=6,units="in",dpi=600)

#
# infomap multilayer
#
fileName <- c(system.file("extdata",  package = "multiweb"))
dn <- list.files(fileName,pattern = "^Kefi2015.*\\.txt$")
g <- readNetwork(dn,fileName, skipColumn = 2)
class(g)
names(g) <- c("Negative","Positive","Trophic")
run_infomap_multi(g,names(g))
plot_troph_level_ggplot(g[[1]], arrow_size = 0.05, label_size = 2)
plot_troph_level_ggplot(g[[2]], arrow_size = 0.05, label_size = 2)
plot_troph_level_ggplot(g[[3]], arrow_size = 0.05, label_size = 2,shorten_factor=0.001)

plot_troph_level(g[[1]])
convert_to_intra_format(g)
names(g)

(gs <- convert_to_supra_adjacency(g,interlayer_weight = 0.4, layer_names= names(g),use_names = TRUE))
ig <- igraph::graph_from_adjacency_matrix(gs$supra_matrix, mode = "directed", weighted = TRUE)
mo <- run_infomap(ig, output_dir = ".")
membership(mo)
mo$codelength
plot_troph_level_ggplot(ig,modules=TRUE,community_obj=mo)

# Result from multilayer Infomap

res_multi <- run_infomap_multi(g, layer_names = names(g),multilayer_relax_rate = 0.4)$communities
multi_modules <- paste0(res_multi$layer, "_", gsub("\\s+", "_", tolower(res_multi$node)))
multi_membership <- setNames(res_multi$module, multi_modules)

# Result from monolayer projection
mono_membership <- membership(mo)
# Clean names in mono_membership for fair comparison
names(mono_membership) <- gsub("\\s+", "_", names(mono_membership))
common_nodes <- intersect(names(multi_membership), names(mono_membership))

comparison <- data.frame(
  Node = common_nodes,
  Multi = multi_membership[common_nodes],
  Mono = mono_membership[common_nodes]
)

comparison$Same <- comparison$Multi == comparison$Mono
print(comparison)
table(comparison$Same)

# To quantitatively compare clustering similarity, use mclust::adjustedRandIndex:

if (!requireNamespace("mclust", quietly = TRUE)) install.packages("mclust")
library(mclust)

adjustedRandIndex(multi_membership[common_nodes], mono_membership[common_nodes])

#
# Centrality
#
calc_centrality(ig,centrality_func = igraph::page_rank)

gs <- convert_to_supra_adjacency(g,interlayer_weight = 0.4, layer_names= names(g),use_names = TRUE,interlayer=FALSE)
ig <- igraph::graph_from_adjacency_matrix(gs$supra_matrix, mode = "directed", weighted = TRUE)
calc_centrality(ig,centrality_func = igraph::page_rank)

calc_svd_entropy_importance(ig)

#
# Shuffling and infoMap
#
names(netData)
result <- generate_shuffled_seq(netData[[29]], modularity_func = run_infomap, shuffle_func = shuffle_network_deg)
# # Extract modularity metrics
# metrics_df <- result$Metrics
#
# # Add a column identifying the original vs shuffled
# metrics_df$Type <- ifelse(metrics_df$Step == 1, "Original", "Shuffled")
#
# # Plot modularity over shuffle steps
# ggplot(metrics_df, aes(x = Step, y = Modularity, color = Type)) +
#   geom_line() +
#   geom_point(size = 2) +
#   labs(title = "Modularity over Network Shuffling",
#        x = "Shuffle Step",
#        y = "Modularity",
#        color = "Network Type") +
#   theme_bw()
#
# ggsave("~/Downloads/degree_shuffling_plot.png", width=8,height=6,units="in",dpi=600)


generate_shuffled_seq(netData[[29]], modularity_func = run_infomap)
generate_shuffled_seq(netData[[23]], modularity_func = run_infomap, shuffle_func = shuffle_network_deg)
generate_shuffled_seq(netData[[23]], modularity_func = run_infomap)

generate_shuffled_seq(netData[[1]], modularity_func = run_infomap, shuffle_func = shuffle_network_ws)
generate_shuffled_seq(netData[[1]], modularity_func = run_infomap)

shuffle_network_deg(netData[[29]],weighted=FALSE)
run_infomap(netData[[29]], output_dir = ".")

require(igraph)
require(NetIndices)
require(ggplot2)
g <-   graph_from_literal( 1 -+ 4 -+ 7,2 -+ 5 -+7, 3-+6-+7, 7-+7, 4+-3, 2-+7, simplify = FALSE)
g1 <-  graph_from_literal( 5 + 4 + 6,3 + 1 + 2, simplify = FALSE)
calc_topological_indices(g)
generate_shuffled_seq(g, modularity_func = run_infomap, shuffle_func = shuffle_network_deg)

plot_troph_level_ggplot(g,modules = TRUE) + ggtitle("Example title")

gs <- convert_to_supra_adjacency(list(g,g1),interlayer_weight = 0.4, layer_names= c("Trophic","Competitive"),use_names = TRUE)
ig <- igraph::graph_from_adjacency_matrix(gs$supra_matrix, mode = "directed", weighted = TRUE)
calc_centrality(ig,centrality_func = igraph::page_rank)
co <- run_infomap(ig)
plot_troph_level_ggplot(ig,modules=TRUE,community_obj=co)

m_list <- list(g,g1)
names(m_list) <- c("Trophic","Competitive")
res <- run_infomap_multi(m_list, output_dir = ".", multilayer_relax_rate = 1 )
plots <- plot_multiplex_modules(m_list, res$communities,y_by_trophic=TRUE, show_labels = TRUE)
cowplot::plot_grid(plotlist = plots, ncol = 2)


# test plot

g <- netData$"Potter Cove"  # example from multiweb

# Run Infomap to get flow per node
infomap_result <- run_infomap(g, return_df = TRUE)
communities_df <- infomap_result$communities

# Match flow to node order
# Make sure V(g)$name matches 'node' in communities_df
# If needed, assign names:
if (is.null(V(g)$name)) {
  V(g)$name <- as.character(1:vcount(g))
}

# Join flow vector: reorder to match igraph node order
node_flow <- communities_df$flow[match(V(g)$name, communities_df$node)]

icom <- run_infomap(g)

# Now plot, using node_weights = Infomap flow
(p <- plot_troph_level_ggraph(
  g,
  modules = TRUE,
  #community_obj = icom,
  #node_weights = node_flow*15,
  vertexSizeFactor = 3,  # increase for visibility
  vertexSizeMin = 1,
  arrow_size = 0.1,
  use_numbers = TRUE
))

(p <- plot_troph_level_ggplot(
  g,
  modules = TRUE,
  community_obj = icom,
  node_weights = NULL,
  vertexSizeFactor = 3,  # increase for visibility
  vertexSizeMin = 1,
  arrow_size = 0.1,
  use_numbers = TRUE
))
plot_troph_level_ggplot(g)
plot_troph_level_ggraph(g)

V(g)$weight <- node_flow*10
(p <- plot_troph_level_ggraph(
  g,
  modules = TRUE,
  community_obj = icom,
  node_weights = NULL,
  vertexSizeFactor = 3,  # increase for visibility
  vertexSizeMin = 1,
  arrow_size = 0.1,
  use_numbers = TRUE
))

pc_i <- calc_interaction_intensity(PotterCove_bm,r_bodymass,r_density,c_bodymass,D)
g <- graph_from_data_frame(pc_i %>% select(resource,consumer, qRC) %>% rename(weight=qRC), directed = TRUE)
(p <- plot_troph_level_ggraph(
  g,
  modules = TRUE,
  weights = NULL,
  vertexSizeFactor = 3,  # increase for visibility
  vertexSizeMin = 1,
  arrow_size = 0.1,
  use_numbers = TRUE
))
coms <- cluster_spinglass(g)
tr <- calc_topological_roles(g,community = coms)
classify_topological_roles(tr,g)

comi <- run_infomap(g, two_level = FALSE)
(p <- plot_troph_level_ggraph(
  g,
  modules = TRUE,
  community_obj = comi,
  vertexSizeFactor = 3,  # increase for visibility
  vertexSizeMin = 1,
  arrow_size = 0.1,
  use_numbers = TRUE
))
tr <- calc_topological_roles(g,community = comi)
classify_topological_roles(tr,g)

comi <- cluster_walktrap(g)

#
# Test calc_stability_threshold
#
g <- generate_niche(40, 0.1)
result <- calc_stability_threshold(g, nsim = 500)
plot_stability_curve(result)
g <- netData[[2]]
calc_topological_indices(g)
result <- calc_stability_threshold(g, nsim = 500)
plot_stability_curve(result)
mean(result$qss_raw_interp[result$qss_raw_interp < 0])
