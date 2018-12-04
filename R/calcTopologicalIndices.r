


#' Calculate topological indices for ecological networks
#'
#' @param ig vector of igraph objects
#' @param ncores number of cores used to compute in parallel, if NULL only one core is used.
#'
#' @return a data.frame with the following fields:
#'
#'  \item{Size:}{Number of species}
#'  \item{Top:}{Number of top predator species}
#'  \item{Basal:}{Number of basal especies}
#'  \item{Omnivory:}{ Number of omnivory species, is the number }
#'  \item{Links:}{ number of interactions}
#'  \item{LD:}{ linkage density}
#'  \item{Connectance:}{ Connectance}
#'  \item{PathLength:}{ average path length}
#'  \item{Clustering:}{ clustering coeficient}
#'  \item{Cannib:}{ number of cannibalistic species}
#'  \item{TLmean:}{ mean trophic level}
#'  \item{Cannib:}{ maximum trophic level}
#'
#'
#' @export
#' @importFrom NetIndices TrophInd
#' @import igraph
#' @importFrom foreach foreach %dopar%
#'
#' @examples
#'
#' topologicalIndices(netData)
#'
#' # Generate a test network
#'
#' g <-   graph_from_literal( 1 -+ 4 -+ 7,2 -+ 5 -+7, 3-+6-+7, 7-+7, 4+-3, simplify = FALSE)
#'
#' topologicalIndices(g)

calcTopologicalIndices <- function(ig,ncores=NULL){

  if(inherits(ig,"igraph")) {
    ig <- list(ig)
  } else if(class(ig[[1]])!="igraph") {
    stop("parameter ig must be an igraph object")
  }

  if(!is.null(ncores)) {
    cn <-parallel::detectCores()
    if(cn>ncores)
      cn <- ncores
    else
      cn <- cn-1
    # cl <- makeCluster(cn,outfile="foreach.log") # Logfile to debug
    cl <- parallel::makeCluster(cn)
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl))
  } else {
    foreach::registerDoSEQ()
  }

  df <-  foreach(g=ig,.combine='rbind',.inorder=FALSE,.packages=c('igraph','NetIndices')) %dopar%
    {

    deg <- degree(simplify(g), mode="out") # calculate the out-degree: the number of predators

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

    pathLength <- average.path.length(g)   # Average path length

    clusCoef <- transitivity(g, type = "global")

    cannib <- sum(which_loop(g))

    tl <- TrophInd(get.adjacency(g,sparse=F))  # Calculate the trophic level

    omn <- sum(round(tl$OI,5)>0)/size        # Omnivory

    data.frame(Size=size,Top=nTop,Basal=nBasal,Omnivory=omn,Links=links, LD=linkDen,Connectance=conn,PathLength=pathLength,
               Clustering=clusCoef, Cannib=cannib, TLmean=mean(tl$TL),TLmax=max(tl$TL))
  }

}


#' Calculate the incoherence index of a food web
#'
#' The incoherence index is based in how the species fit in discrete trophic levels
#' when Q is closer to 0 more coherent and stable is a food web, and the less omnivory it has.
#' It calculates the trophic level using the package NetIndices, or optionally the trophic levels could be passed
#' by parameter as a vector but in that case that the same vector will be used for all the networks.
#'
#' Based on:
#'
#' @references Johnson, S., Domínguez-García, V., Donetti, L., & Muñoz, M. A. (2014). Trophic coherence determines food-web stability. Proceedings of the National Academy of Sciences , 111(50), 17923–17928. https://doi.org/10.1073/pnas.1409077111
#'
#' @references Johnson, S., & Jones, N. S. (2017). Looplessness in networks is linked to trophic coherence. Proceedings of the National Academy of Sciences, 114(22), 5618–5623. https://doi.org/10.1073/pnas.1613786114
#'
#'
#' @param g an igraph object
#' @param ti trophic level vector if not supplied is calculated internally
#' @param ncores number of cores used to compute in parallel, if NULL one core is used.
#'
#' @return a data.frame with the following fields
#'
#'  \item{Q}{incoherence (0=coherent)}
#'  \item{rQ}{ratio of Q with expected Q under null expectation of a random network given N=nodes L=links B=basal nodes Lb=basal links}
#'  \item{mTI}{mean trophic level}
#'  \item{rTI}{ratio of mTI with expected TI under the same null model expectation than Q}
#
#'
#' @export
#'
#' @examples
#'
#' calcIncoherence(netData[[1]])
#'
#'
#' @importFrom NetIndices TrophInd
#' @importFrom igraph     V degree get.adjacency vcount ecount
calcIncoherence <- function(ig,ti=NULL,ncores=NULL) {

  if(inherits(ig,"igraph")) {
    ig <- list(ig)
  } else if(class(ig[[1]])!="igraph") {
    stop("parameter ig must be an igraph object")
  }

  if(!is.null(ncores)) {
    cn <-parallel::detectCores()
    if(cn>ncores)
      cn <- ncores
    else
      cn <- cn-1

    # cl <- makeCluster(cn,outfile="foreach.log") # Logfile to debug
    cl <- parallel::makeCluster(cn)
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl))
  } else {
    foreach::registerDoSEQ()
  }

  df <-  foreach(g=ig,.combine='rbind',.inorder=FALSE,.packages=c('igraph','NetIndices')) %dopar%
  {

    if(is.null(ti) )
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
}


#' Calculation of Modularity and Small-world-ness z-scores
#'
#' The function calculates modularity, number of groups and small-world-ness z-scores and 99\% CI intervals
#' using as null model Erdos-Renyi random networks with the same number of links and species.
#'
#'
#' @references
#' Marina, T. I., Saravia, L. A., Cordone, G., Salinas, V., Doyle, S. R., & Momo, F. R. (2018). Architecture of marine food webs: To be or not be a ‘small-world.’ PLoS ONE, 13(5), 1–13. https://doi.org/10.1371/journal.pone.0198217
#'
#' @param g  igraph object
#' @param nullDist list of igraph object with the null model simulations
#' @param sLevel significance level to calculate CI (two tails)
#' @param ncores number of cores to use paralell computation, if NULL no parallel computation is used
#'
#'
#' @return a data frame with indices z-scores and CI
#'
#'  \item{Clustering}{ Clustering coefficient, measures the average fraction of pairs of neighbors of a node that are also neighbors of each other}
#'  \item{PathLength}{ Mean of the shortest paths between all pair of vertices }
#'  \item{Modularity}{ modularity measures how separated are different groups from each other, the algorithm \code{cluster_spinglass} was used to obtain the groups}
#'  \item{zCC,zCP,zMO}{Z-scores of Clustering,PathLength and Modularity with respect to a random Erdos-Renyi null model}
#'  \item{CClow,CChigh,CPlow,CPhigh,MOlow,MOhigh}{sLevel confidence intervals}
#'  \item{SWness,SWnessCI}{ Small-world-ness and it CI value}
#'  \item{isSW,isSWness}{ Logical variable signalling if the network is Small-world by the method of Marina 2018 or the method of Humprhies & Gurney 2008 }

#' @export
#'
#' @importFrom foreach foreach %dopar%
#' @importFrom igraph transitivity average.path.length cluster_spinglass

#' @examples
#'
#' calcModularitySWnessZScore(netData[[1]])

calcModularitySWnessZScore<- function(g, nullDist,sLevel=0.01,ncores=NULL){

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

  if(!is.null(ncores)) {
    cn <-parallel::detectCores()
    if(cn>ncores)
      cn <- ncores
    else
      cn <- cn-1

    # cl <- makeCluster(cn,outfile="foreach.log") # Logfile to debug
    cl <- parallel::makeCluster(cn)
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl))
  } else {
    foreach::registerDoSEQ()
  }

  ind <- data.frame()

  ind <- foreach(i=1:length(nullDist),.combine='rbind',.inorder=FALSE,.packages='igraph') %dopar%
  {
    m<-cluster_spinglass(nullDist[[i]])
    modl <- m$modularity
    clus.coef <- transitivity(nullDist[[i]], type="Global")
    cha.path  <- average.path.length(nullDist[[i]])
    data.frame(modularity=modl,clus.coef=clus.coef,cha.path=cha.path)
  }

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



  m<-cluster_spinglass(g)

  zmo <- (m$modularity - mean(ind$modularity))/sd(ind$modularity)

  mSW <- mean(t$Clustering/mcc*mcp/t$PathLength)
  mCI <- 1+(qSW[2]-qSW[1])/2

  isLowInsideCPL <- t$PathLength<=qcp[2]
  isGreaterCC <- t$Clustering>qcc[2]
  isSW       <- (isLowInsideCPL & isGreaterCC)
  isSWness   <- (mSW>mCI)

  return(data.frame(Clustering=t$Clustering, PathLength= t$PathLength, Modularity=m$modularity,
    zCC=zcc,zCP=zcp,zMO=zmo,
    CClow=qcc[1],CChigh=qcc[2],CPlow=qcp[1],CPhigh=qcp[2],MOlow=qmo[1],MOhigh=qmo[2],
    SWness=mSW,SWnessCI=mCI,
    isSW=isSW,isSWness=isSWness))
}
