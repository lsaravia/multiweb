


#' Calculate topological indices for ecological networks. The mean trophic level
#' and Omnivory and the Level of omnivory are calculated with the function [NetIndices::TrophInd()]
#'
#' @param ig vector of igraph objects
#' @param ncores number of cores used to compute in parallel, if 0 sequential processing is used.
#'
#' @return a data.frame with the following fields:
#'
#'  * Network: Name of the network if the list have names
#'  * Size: Number of species
#'  * Top: Number of top predator species
#'  * Basal: Number of basal especies
#'  * Omnivory:  Proportion of omnivorous species
#'  * Links:  number of interactions
#'  * LD:     linkage density
#'  * Connectance:  directed Connectance
#'  * PathLength:   average path length
#'  * Clustering:  clustering coeficient
#'  * Cannib:  number of cannibalistic species
#'  * TLmean:  mean trophic level
#'  * TLmax:   maximum trophic level
#'  * LOmnivory:  Level of omnivory, quantiﬁes mean of the variety in trophic levels of the preys of a consumer
#'  * Components:  number of weakly connected components
#'  * Vulnerability:  mean of number of consumers per prey
#'  * VulSD:  the standard deviation of normalized Vulnerability
#'  * Generality:  mean number of prey per consumer
#'  * GenSD:  the standard deviation of normalized Generality
#'
#'
#' @export
#' @importFrom NetIndices TrophInd
#' @import igraph
#' @importFrom foreach foreach %dopar%
#' @importFrom doFuture registerDoFuture
#' @importFrom future sequential multisession
#' @importFrom dplyr mutate %>% select
#'
#' @aliases calcTopologicalIndices
#'
#' @examples
#'
#' calc_topological_indices(netData)
#'
#' # Generate a test network
#'
#' g <-   graph_from_literal( 1 -+ 4 -+ 7,2 -+ 5 -+7, 3-+6-+7, 7-+7, 4+-3, simplify = FALSE)
#'
#' calc_topological_indices(g)
#'
calc_topological_indices <- function(ig,ncores=0){

  if(inherits(ig,"igraph")) {
    ig <- list(ig)
  } else if(class(ig[[1]])!="igraph") {
    stop("parameter ig must be an igraph object")
  }

  registerDoFuture()
  if(ncores) {
    cn <- future::availableCores()
    if(ncores>cn)
      ncores <- cn
    future::plan(multisession, workers=ncores)
    on.exit(future::plan(sequential))
  } else {
    future::plan(sequential)
  }


  df <-  foreach(g=ig,.combine='rbind',.inorder=FALSE,.packages=c('igraph','NetIndices')) %dopar%
    {

    deg <- degree(igraph::simplify(g), mode="out") # calculate the out-degree: the number of predators

    V(g)$outdegree <-  deg

    nTop <- length(V(g)[outdegree==0]) # Top predators do not have predators

    deg <- degree(g, mode="in") # calculate the in-degree: the number of preys

    V(g)$indegree <-  deg

    nBasal <- length(V(g)[indegree==0]) # Basal species do not have preys

    #vcount(g)-nTop-nBasal

    size <- vcount(g)

    links <- ecount(g)

    linkDen <- links/size          # Linkage density

    conn <- links/size^2           # Connectance

    pathLength <- mean_distance(g)   # Average path length

    clusCoef <- transitivity(g, type = "global")

    cannib <- sum(which_loop(g))

    a <- as_adjacency_matrix(g,sparse=F)

    tl <- TrophInd(a)  # Calculate the trophic level

    omn <- sum(round(tl$OI,5)>0)/size        # Omnivory proportion

    lomn <- mean(tl$OI)                      # Level of omnivory

    vulnerability <- links / (size - nTop)    # l/(nb + ni)
    generality    <- links / (size - nBasal)  # l/(nt + ni)
    vulSD         <- sd(apply(a, 1, sum)*1/linkDen)
    genSD         <- sd(apply(a, 2, sum)*1/linkDen)

    data.frame(Size=size,Top=nTop,Basal=nBasal,Omnivory=omn,Links=links, LD=linkDen,Connectance=conn,PathLength=pathLength,
               Clustering=clusCoef, Cannib=cannib, TLmean=mean(tl$TL),TLmax=max(tl$TL),LOmnivory=lomn,Components=components(g)$no,
               Vulnerability=vulnerability,VulSD=vulSD,Generality=generality,GenSD=genSD)
    }
  if( !is.null(names(ig)) )
    df <- df %>% mutate(Network = names(ig)) %>% dplyr::select(Network, everything())
  return(df)
}

#' @export
calcTopologicalIndices <- function(ig,ncores=0){
    calc_topological_indices(ig,ncores)
}


#' Calculate the incoherence index of a food web
#'
#' The incoherence index is based in how the species fit in discrete trophic levels
#' when Q is closer to 0 more coherent and stable is a food web, and the less omnivory it has.
#' It calculates the trophic level using the package NetIndices.
#'
#' Based on:
#'
#' @references Johnson, S., Domínguez-García, V., Donetti, L., & Muñoz, M. A. (2014). Trophic coherence determines food-web stability. Proceedings of the National Academy of Sciences , 111(50), 17923–17928. https://doi.org/10.1073/pnas.1409077111
#'
#' @references Johnson, S., & Jones, N. S. (2017). Looplessness in networks is linked to trophic coherence. Proceedings of the National Academy of Sciences, 114(22), 5618–5623. https://doi.org/10.1073/pnas.1613786114
#'
#'
#' @param ig an igraph object or a list of igraph objects
#' @param ncores number of cores used to compute in parallel, if 0 sequential processing is used.
#'
#' @return a data.frame with the following fields
#'
#'  \item{Q}{incoherence (0=coherent)}
#'  \item{rQ}{ratio of Q with expected Q under null expectation of a random network given N=nodes L=links B=basal nodes Lb=basal links}
#'  \item{mTI}{mean trophic level}
#'  \item{rTI}{ratio of mTI with expected TI under the same null model expectation than Q}
#
#' @aliases calcIncoherence
#'
#' @export
#'
#' @examples
#'
#' calc_incoherence(netData[[1]])
#'
#'
#' @importFrom NetIndices TrophInd
#' @importFrom igraph     V degree get.adjacency vcount ecount
#' @importFrom foreach foreach %dopar%
#' @importFrom doFuture registerDoFuture
#' @importFrom future sequential multisession

calc_incoherence <- function(ig,ncores=0) {

  if(inherits(ig,"igraph")) {
    ig <- list(ig)
  } else if(class(ig[[1]])!="igraph") {
    stop("parameter ig must be an igraph object")
  }

  if(ncores) {
    cn <- future::availableCores()
    if(ncores>cn)
      ncores <- cn
    future::plan(multisession, workers=ncores)
    on.exit(future::plan(sequential))
  } else {
    future::plan(sequential)
  }


  # if(!is.null(ncores)) {
  #   cn <-parallel::detectCores()
  #   if(cn>ncores)
  #     cn <- ncores
  #   else
  #     cn <- cn-1
  #
  #   # cl <- makeCluster(cn,outfile="foreach.log") # Logfile to debug
  #   cl <- parallel::makeCluster(cn)
  #   doParallel::registerDoParallel(cl)
  #   on.exit(parallel::stopCluster(cl))
  # } else {
  #   foreach::registerDoSEQ()
  # }

  df <-  foreach(g=ig,.combine='rbind',.inorder=FALSE,.packages=c('igraph','NetIndices')) %dopar%
  {
    ti<-TrophInd(get.adjacency(g,sparse=FALSE))
    v <- ti$TL
    z <- round(outer(v,v,'-'),8);
    A <- get.adjacency(g,sparse = FALSE)
    xx <- A>0
    x <- (A*t(z))[xx]
    Q <- round(sqrt(sum((x-1)^2)/ecount(g) ),8)

    basal <- which(round(v,8)==1)
    bedges <- sum(degree(g,basal,mode='out'))
    mTI <- mean(v)
    mK <- mean(degree(g,mode='out'))
    eTI <- 1+(1-length(basal)/vcount(g))*ecount(g)/bedges
    eQ <- sqrt(ecount(g)/bedges-1)
    data.frame(Q=Q,rQ=Q/eQ,mTI=mTI,rTI=mTI/eTI)
  }
  return(df)
}

#' @export
calcIncoherence <- function(ig,ncores=0){
  calc_incoherence(ig,ncores)
  }

#' Calculation of Modularity and Small-world-ness z-scores
#'
#' The function calculates modularity, number of groups and small-world-ness z-scores and 99\% CI intervals
#' using as null model the list of networks in the nullDist parameters. Modularity is calculated using the [igraph::cluster_spinglass()]
#' if the parameter weights is NULL the atribute "weigths" is used, or if it has the name of an network attribute, that is used as a weigth
#' to build the modules, when this parameter is NA then no weigth is used.
#' Only works for one component networks
#'
#' @references
#' Marina, T. I., Saravia, L. A., Cordone, G., Salinas, V., Doyle, S. R., & Momo, F. R. (2018). Architecture of marine food webs: To be or not be a ‘small-world.’ PLoS ONE, 13(5), 1–13. https://doi.org/10.1371/journal.pone.0198217
#'
#' @param g  igraph object
#' @param nullDist list of igraph object with the null model simulations
#' @param sLevel significance level to calculate CI (two tails)
#' @param ncores number of cores to use paralell computation, if 0 sequential processing is used.
#' @param weights The weights of the edges. Either a numeric vector or NULL or NA. If it is null and the input graph has a ‘weight’ edge attribute
#'                then that will be used. If NULL and no such attribute is present then the edges will have equal weights.
#'                Set this to NA if the graph was a ‘weight’ edge attribute, but you don't want to use it for community detection.
#'
#'
#' @return a list with two data frames: one with indices z-scores and CI
#'
#'  \item{Clustering}{ Clustering coefficient, measures the average fraction of pairs of neighbors of a node that are also neighbors of each other}
#'  \item{PathLength}{ Mean of the shortest paths between all pair of vertices }
#'  \item{Modularity}{ modularity measures how separated are different groups from each other, the algorithm \code{cluster_spinglass} was used to obtain the groups}
#'  \item{zCC,zCP,zMO}{Z-scores of Clustering,PathLength and Modularity with respect to a random Erdos-Renyi null model}
#'  \item{CClow,CChigh,CPlow,CPhigh,MOlow,MOhigh}{sLevel confidence intervals}
#'  \item{SWness,SWnessCI}{ Small-world-ness and it CI value}
#'  \item{isSW,isSWness}{ Logical variable signalling if the network is Small-world by the method of Marina 2018 or the method of Humprhies & Gurney 2008 }
#'
#'  Another data.frame with the values calculated for the nullDist.
#'
#' @export
#'
#' @importFrom igraph transitivity average.path.length cluster_spinglass
#' @importFrom future.apply future_lapply
#' @importFrom future sequential multisession
#'
#' @examples
#' \dontrun{
#' nullg <- generateERbasal(netData[[1]],10)
#' calcModularitySWnessZScore(netData[[1]],nullg)
#' }
#'
calc_modularity_swness_zscore<- function(g, nullDist,sLevel=0.01,ncores=0,weights=NA){

  if(!is_igraph(g))
    stop("Parameter g must be an igraph object")

  t <- calcTopologicalIndices(g)

  if(length(nullDist)<5)
    stop("nullDist: There has to be more than 5 elements in the list")

  if(any(sapply(nullDist, function(g) components(g)$no>1)))
    stop("nullDist: one or more igraph object from nullDist have more than one component")

  # nullDist <- lapply(1:nsim, function (x) {
  #   e <- sample_gnm(t$Size, t$Links, directed = TRUE)
  #
  #   # Check that the ER networks has only one connected component
  #   #
  #   while(components(e)$no>1)
  #     e <- erdos.renyi.game(t$Size, t$Links, type="gnm",directed = TRUE)
  #
  #   return(e) }
  # )

  if(ncores) {
    cn <- future::availableCores()
    if(ncores>cn)
      ncores <- cn
    future::plan(multisession, workers=ncores)
    on.exit(future::plan(sequential))
  } else {
    future::plan(sequential)
  }

  ind <- data.frame()

  ind <- future_lapply(1:length(nullDist), function(i)
  {
    m<-cluster_spinglass(nullDist[[i]],weights=weights)
    modl <- m$modularity
    clus.coef <- transitivity(nullDist[[i]], type="Global")
    cha.path  <- average.path.length(nullDist[[i]])
    data.frame(modularity=modl,clus.coef=clus.coef,cha.path=cha.path)
  }, future.seed = TRUE)
  ind <- bind_rows(ind)

  ind$gamma <- t$Clustering/ind$clus.coef
  ind$lambda <- t$PathLength/ind$cha.path
  ind$SWness <- ind$gamma/ind$lambda

  # 99% confidence interval
  #
  sLevel <- sLevel/2
  qSW <- quantile(ind$SWness,c(sLevel,1-sLevel),na.rm = TRUE)
  qmo <- quantile(ind$modularity,c(sLevel,1-sLevel))

  mcc <- mean(ind$clus.coef)
  mcp <- mean(ind$cha.path)

  zcc <- (t$Clustering-mcc)/sd(ind$clus.coef)
  zcp <- (t$PathLength-mcp)/sd(ind$cha.path)
  qcc <- quantile(ind$clus.coef,c(sLevel,1-sLevel),na.rm = TRUE)
  qcp <- quantile(ind$cha.path,c(sLevel,1-sLevel),na.rm = TRUE)



  m<-cluster_spinglass(g, weights=weights)

  zmo <- (m$modularity - mean(ind$modularity))/sd(ind$modularity)

  mSW <- mean(t$Clustering/mcc*mcp/t$PathLength)
  mCI <- 1+(qSW[2]-qSW[1])/2

  isLowInsideCPL <- t$PathLength<=qcp[2]
  isGreaterCC <- t$Clustering>qcc[2]
  isSW       <- (isLowInsideCPL & isGreaterCC)
  isSWness   <- (mSW>mCI)

  return(list(da=data.frame(Clustering=t$Clustering, PathLength= t$PathLength, Modularity=m$modularity,
    zCC=zcc,zCP=zcp,zMO=zmo,
    CClow=qcc[1],CChigh=qcc[2],CPlow=qcp[1],CPhigh=qcp[2],MOlow=qmo[1],MOhigh=qmo[2],
    SWness=mSW,SWnessCI=mCI,
    isSW=isSW,isSWness=isSWness),sims=ind))
}

#' @export
calcModularitySWnessZScore<- function(g, nullDist,sLevel=0.01,ncores=0){
  calc_modularity_swness_zscore(g, nullDist,sLevel,ncores)
}


#' Calculation Small-world-ness z-scores
#'
#' The function calculates small-world-ness z-scores and 99\% CI intervals,
#' using as null model the list of networks in the nullDist parameter.
#'
#' @references
#' Marina, T. I., Saravia, L. A., Cordone, G., Salinas, V., Doyle, S. R., & Momo, F. R. (2018). Architecture of marine food webs: To be or not be a ‘small-world.’ PLoS ONE, 13(5), 1–13. https://doi.org/10.1371/journal.pone.0198217
#'
#' @param g  igraph object
#' @param nullDist list of igraph object with the null model simulations
#' @param sLevel significance level to calculate CI (two tails)
#' @param ncores number of cores to use paralell computation, if 0 sequential processing is used.
#'
#'
#' @return a list with two data frames: one with indices z-scores and CI
#'
#'  \item{Clustering}{ Clustering coefficient, measures the average fraction of pairs of neighbors of a node that are also neighbors of each other}
#'  \item{PathLength}{ Mean of the shortest paths between all pair of vertices }
#'  \item{Modularity}{ modularity measures how separated are different groups from each other, the algorithm \code{cluster_spinglass} was used to obtain the groups}
#'  \item{zCC,zCP,zMO}{Z-scores of Clustering,PathLength and Modularity with respect to `nullDist` models}
#'  \item{CClow,CChigh,CPlow,CPhigh,MOlow,MOhigh}{sLevel confidence intervals}
#'  \item{SWness,SWnessCI}{ Small-world-ness and it CI value}
#'  \item{isSW,isSWness}{ Logical variable signalling if the network is Small-world by the method of Marina 2018 or the method of Humprhies & Gurney 2008 }
#'
#'  Another data.frame with the values calculated for the nullDist.
#'
#' @export
#'
#' @importFrom igraph transitivity average.path.length
#' @importFrom foreach foreach %dopar%
#' @importFrom doFuture registerDoFuture
#' @importFrom future sequential multisession
#'
#' @examples
#' \dontrun{
#' nullg <- generateERbasal(netData[[1]],10)
#' calc_swness_zscore(netData[[1]],nullg)
#' }
#'
calc_swness_zscore<- function(g, nullDist,sLevel=0.01,ncores=0,weights=NA){

  if(!is_igraph(g))
    stop("Parameter g must be an igraph object")

  t <- calcTopologicalIndices(g)

  if(length(nullDist)<5)
    stop("nullDist: There has to be more than 5 elements in the list")

  registerDoFuture()
  if(ncores) {
    cn <- future::availableCores()
    if(ncores>cn)
      ncores <- cn
    future::plan(multisession, workers=ncores)
    on.exit(future::plan(sequential))
  } else {
    future::plan(sequential)
  }

  ind <- data.frame()

  ind <- foreach(i=1:length(nullDist),.combine='rbind',.inorder=FALSE,.packages='igraph') %dopar%
    {
      clus.coef <- transitivity(nullDist[[i]], type="Global")
      cha.path  <- average.path.length(nullDist[[i]])
      data.frame(clus.coef=clus.coef,cha.path=cha.path)
    }

  ind$gamma <- t$Clustering/ind$clus.coef
  ind$lambda <- t$PathLength/ind$cha.path
  ind$SWness <- ind$gamma/ind$lambda

  # 99% confidence interval
  #
  sLevel <- sLevel/2
  qSW <- quantile(ind$SWness,c(sLevel,1-sLevel),na.rm = TRUE)

  mcc <- mean(ind$clus.coef)
  mcp <- mean(ind$cha.path)

  zcc <- (t$Clustering-mcc)/sd(ind$clus.coef)
  zcp <- (t$PathLength-mcp)/sd(ind$cha.path)
  qcc <- quantile(ind$clus.coef,c(sLevel,1-sLevel),na.rm = TRUE)
  qcp <- quantile(ind$cha.path,c(sLevel,1-sLevel),na.rm = TRUE)

  mSW <- mean(t$Clustering/mcc*mcp/t$PathLength)
  mCI <- 1+(qSW[2]-qSW[1])/2

  isLowInsideCPL <- t$PathLength<=qcp[2]
  isGreaterCC <- t$Clustering>qcc[2]
  isSW       <- (isLowInsideCPL & isGreaterCC)
  isSWness   <- (mSW>mCI)

  return(list(da=data.frame(Clustering=t$Clustering, PathLength= t$PathLength,
                            zCC=zcc,zCP=zcp,
                            CClow=qcc[1],CChigh=qcc[2],CPlow=qcp[1],CPhigh=qcp[2],
                            SWness=mSW,SWnessCI=mCI,
                            isSW=isSW,isSWness=isSWness),sims=ind))
}




#' Calculate Topological Roles Within Network Modules
#'
#' This function identifies the topological roles of nodes (e.g., species) in an ecological or interaction network based on their position within and among modules.
#' Modules (or communities) are detected using the \code{\link[igraph]{cluster_spinglass}} algorithm by default, unless an alternative community object is provided.
#'
#' Topological roles are quantified using two parameters:
#' \itemize{
#'   \item \eqn{z} (within-module degree z-score): quantifies how well-connected a node is to others within its own module, relative to other members of the same module.
#'   \item \eqn{PC} (participation coefficient): quantifies how evenly a node’s links are distributed across different modules.
#' }
#'
#' When using the Spinglass algorithm, which is stochastic, multiple runs can be performed by setting the \code{nsim} parameter to repeat the procedure.
#'
#' @param g An \code{igraph} object representing the network.
#' @param nsim Integer. Number of independent runs of the Spinglass community detection algorithm (default is 1).
#' @param ncores Integer. Number of cores to use for parallel computation. Set to 0 for sequential processing.
#' @param community Optional. A community detection object (class \code{communities}) from an alternative algorithm. If \code{NULL}, \code{\link[igraph]{cluster_spinglass}} is used.
#'
#' @return A data frame with the following columns:
#' \describe{
#'   \item{node}{Node ID.}
#'   \item{within_module_degree}{Within-module degree z-score (\eqn{z}).}
#'   \item{among_module_conn}{Participation coefficient (\eqn{PC}).}
#' }
#'
#' @references
#' \enumerate{
#'   \item Guimerà, R. & Nunes Amaral, L.A. (2005). Functional cartography of complex metabolic networks. \emph{Nature}, 433, 895–900.
#'   \item Kortsch, S., Primicerio, R., Fossheim, M., Dolgov, A.V., & Aschan, M. (2015). Climate change alters the structure of Arctic marine food webs due to poleward shifts of boreal generalists. \emph{Proceedings of the Royal Society B: Biological Sciences}, 282(1814), 20151546. \doi{10.1098/rspb.2015.1546}
#' }
#'
#' @export
#'
#' @import igraph
#' @importFrom foreach foreach %dopar%
#' @importFrom future.apply future_lapply
#' @importFrom future sequential multisession
#'
#' @examples
#' #' \dontrun{
#'
#' g <- netData[[2]]
#'
#' tp <- calc_topological_roles(g,nsim=10)
#'
#' # using a community object
#'
#' m <- cluster_walktrap(g)
#'
#' tp <- calc_topological_roles(g,community=m)
#' }

calc_topological_roles <- function(g,nsim=1,ncores=0,community=NULL)
{
  if(!is_igraph(g))
    stop("Parameter g must be an igraph object")

  # Check if community object is provided and its class
  if (!is.null(community)){
    if (class(community) != "communities") {
      stop("Provided community object must be of class 'communities'")
    } else {
      nsim <- 1
      ncores <- 0
    }
  }
  toRol <- data.frame()


  if(ncores) {
    cn <- future::availableCores()
    if(ncores>cn)
      ncores <- cn
    future::plan(multisession, workers=ncores)
    on.exit(future::plan(sequential))
  } else {
    future::plan(sequential)
  }
  toRol <- do.call(rbind, future_lapply(seq_len(nsim), function(idx){

    # within-module degree
    #
    # Standarized Within module degree z-score
    #
    if (is.null(community)) {
      m <- cluster_spinglass(g, weights = NA)
      spingB.mem <- m$membership
    } else {
      spingB.mem <- membership(community)
    }

    l<-vector()
    memMod<-vector()

    for (i in 1:vcount(g)){

      sp.in.deg <- V(g)[.nei(i, "in")]
      sp.out.deg<- V(g)[.nei(i, "out")]
      mem.sp<-spingB.mem[i]
      k<- length(which(spingB.mem[c(sp.in.deg, sp.out.deg)]==mem.sp))
      mem<- which(spingB.mem==mem.sp)

      for (m in 1:length(mem)){
        mem.in.deg <- V(g)[.nei(mem[m], "in")]
        mem.out.deg<- V(g)[.nei(mem[m], "out")]
        memMod.id<- length(which(spingB.mem[c(mem.in.deg, mem.out.deg)]==mem.sp))
        memMod[m]<- memMod.id
      }

      k.ave<- mean(memMod)
      k.sd<- sd(memMod)
      l[i]<- ifelse(k.sd == 0, 0, (k-k.ave)/k.sd)
    }

    # among module connectivity
    r<- vector()
    for (i in 1:vcount(g)){

      d<-degree(g)[i]
      sp.in.deg <- V(g)[.nei(i, "in")]
      sp.out.deg<- V(g)[.nei(i, "out")]
      mod.sp<-table(spingB.mem[c(sp.in.deg, sp.out.deg)])
      mod.no<-as.numeric(names(mod.sp))
      mod<-rep(0, length(unique(spingB.mem)))
      mod[mod.no]<-c(mod.sp)
      r[i]<- 1-(sum((mod/d)^2))
    }
    return(data.frame(node=1:vcount(g),within_module_degree=l, among_module_conn=r))

  }, future.seed = TRUE ))

  return(toRol)
}


#' Classify and plot topological roles
#'
#' This function classifies species into four topological roles based on their within-module degree (dz) and among-module connectivity (PC), following the approach of Kortsch (2015). The classification is as follows:
#' - **modhub**: Species with high within-module degree (dz >= 2.5) and low among-module connectivity (PC < 0.625).
#' - **hubcon**: Species with high within-module degree (dz >= 2.5) and high among-module connectivity (PC >= 0.625).
#' - **modspe**: Species with low within-module degree (dz < 2.5) and low among-module connectivity (PC < 0.625).
#' - **modcon**: Species with low within-module degree (dz < 2.5) and high among-module connectivity (PC >= 0.625).
#'
#' The function also generates a ggplot visualization of the species' topological roles, with nodes colored according to their module membership.
#'
#' @references
#' 1. Kortsch, S., Primicerio, R., Fossheim, M., Dolgov, A. V., & Aschan, M. (2015). Climate change alters the structure of arctic marine food webs due to poleward shifts of boreal generalists.
#'    Proceedings of the Royal Society B: Biological Sciences, 282(1814), 20151546. https://doi.org/10.1098/rspb.2015.1546
#' 2. Saravia, L. A., Marina, T. I., Kristensen, N. P., De Troch, M., & Momo, F. R. (2022). Ecological network assembly: How the regional metaweb influences local food webs.
#'    Journal of Animal Ecology, 91(3), 630–642. https://doi.org/10.1111/1365-2656.13652
#'
#' @param tRoles A data frame of calculated topological roles obtained from [calc_topological_roles()].
#' @param g An igraph network object.
#' @param community An igraph community object used to determine module membership for plotting. If NULL, the function calculates community structure using `cluster_spinglass()`.
#'
#' @return A list containing:
#' - `hub_conn`: A data frame with node classifications and topological role values.
#' - `gp`: A ggplot object visualizing the topological roles.
#'
#' @export
#'
#' @import igraph ggplot2
#' @import ggrepel
#' @importFrom dplyr %>% group_by summarise_all inner_join select mutate arrange
#'
#' @examples
#' \dontrun{
#'
#' g <- netData[[2]]
#'
#' tp <- calc_topological_roles(g,nsim=10)
#'
#' classify_topological_roles(tp,g)
#' }
classify_topological_roles <- function(tRoles, g, community = NULL) {

  if (is.null(community)) {
    community <- cluster_spinglass(g, weights = NA)
  }

  spingB.mem <- community$membership
  tRoles <- tRoles %>% group_by(node) %>% summarise_all(mean)

  l <- tRoles$within_module_degree
  r <- tRoles$among_module_conn

  # Data for ggplot
  plot_data <- data.frame(
    node = tRoles$node,
    node_name = vertex_attr(g, "name", index = tRoles$node),  # Include node names
    within_module_degree = l,
    among_module_conn = r,
    module = factor(spingB.mem)
  )

  gp <- ggplot(plot_data, aes(x = among_module_conn, y = within_module_degree, color = module)) +
    geom_point(size = 4, alpha=0.9) +
    geom_vline(xintercept = 0.625, linetype = "dashed", color = "gray") +
    geom_hline(yintercept = 2.5, linetype = "dashed", color = "gray") +
    #scale_color_brewer(palette = "Dark2") +
    scale_color_viridis_d(option = "D") +  # <- 🔑
    labs(x = "Among module connectivity", y = "Within module degree", color = "Module") +
    theme_bw() +
    geom_text_repel(aes(label=node_name))+
    theme(text = element_text(size = 14))


  classify_nodes <- function(l, r) {
    if (l > 2.5 & r <= 0.625) return("modhub")
    if (l > 2.5 & r > 0.625) return("hubcon")
    if (l <= 2.5 & r <= 0.625) return("modspe")
    if (l <= 2.5 & r > 0.625) return("modcon")
    return(NA)
  }

  plot_data$type <- mapply(classify_nodes, plot_data$within_module_degree, plot_data$among_module_conn)

  hub_conn <- plot_data %>% select(node, node_name, type) %>%
    inner_join(tRoles, by = "node") %>%
    mutate(type = factor(type, levels = c("hubcon", "modhub", "modcon", "modspe"))) %>%
    arrange(type)

  return(list(hub_conn=hub_conn,gp=gp))
}


