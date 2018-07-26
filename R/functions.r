# Function definitions



#' Read ecological networks in CSV  or tab separated file format as edge list or adyacency matrix
#'
#' @param fileName vector of fileNames with the networks
#' @param filePath path of the files NULL by default
#' @param fhead TRUE if the files have header fields, FALSE otherwise.
#'
#' @return an igraph object if there is only one file or a list of igraph objects named after the list without extension
#' @export
#'
#' @examples
#'
#' # Reads a network in edge list (interaction list) format, with predators as the first column
#' #
#' fileName <- system.file("extdata", "WeddellSea_FW.csv", package = "EcoNetwork")
#' g <- readNetwork(fileName)
#'
#' # Reads a network in adyacency matrix format, with predators as columns
#' #
#' fileName <- system.file("extdata", "BarentsBoreal_FW.csv", package = "EcoNetwork")
#' g <- readNetwork(fileName)
#'
#' # Read a vector of files
#' #
#'\dontrun{
#' dn <- list.files("inst/extdata",pattern = "^.*\\.csv$")
#' netData <- readNetwork(dn,"inst/extdata")
#'}

readNetwork <- function(fileName,filePath=NULL,fhead=TRUE){

  fn <-  if(!is.null(filePath)) paste0(filePath,"/",fileName) else fileName

  g <- lapply(fn, function(fname){

    fe <- tools::file_ext(fname)
    if(fe=="csv")
      web <- read.csv(fname,  header = fhead,check.names = F,stringsAsFactors = FALSE)
    else {

      web <- read.delim(fname,  header = fhead,check.names = F,stringsAsFactors = FALSE)
      if(ncol(web)==1)
        web <- read.table(fname,header = fhead,check.names = F,stringsAsFactors = FALSE)

    }

    if( ncol(web)==2 ){
      web <- web[,c(2,1)]

      g <- igraph::graph_from_data_frame(web)

    } else {
      if( (ncol(web)-1) == nrow(web)  ) {                   # The adjacency matrix must be square
        g <- igraph::graph_from_adjacency_matrix(as.matrix(web[,2:ncol(web)]))

      } else {
        g <- NULL
        warning("Invalid file format: ",fname)
      }
    }

  })

  names(g) <- tools::file_path_sans_ext(fileName)

  if(length(g)==1)
    g <- g[[1]]

  return(g)
}



#' Plot ecological network organized by trophic level, with node size determined by the node degree
#'
#' @param ig igraph object
#' @param vertexLabel logical plot vertex labels
#' @param vertexSizeFactor numeric factor to determine the size of the label with degree
#' @param tk TRUE generate an interactive plot using tkplot
#'
#' @return a plot
#' @export
#'
#'
#' @importFrom NetIndices TrophInd
#' @importFrom igraph     V degree get.adjacency
#' @importFrom RColorBrewer brewer.pal
#' @examples
#'
#' plotTrophLevel(netData[[1]])
plotTrophLevel <- function(g,vertexLabel=FALSE,vertexSizeFactor=5,tk=FALSE){

  deg <- degree(g, mode="all") # calculate the degree: the number of edges
  # or interactions

  V(g)$size <- log10(deg)*vertexSizeFactor+vertexSizeFactor    # add node degrees to igraph object

  V(g)$frame.color <- "white"    # Specify plot options directly on the object

  V(g)$color <- "orange"         #

  if(!vertexLabel)
    V(g)$label <- NA

  tl <- TrophInd(get.adjacency(g,sparse=F))  # Calculate the trophic level

  # Layout matrix to specify the position of each vertix
  # Rows equal to the number of vertices (species)
  lMat <-matrix(
    nrow=vcount(g),
    ncol=2
  )

  lMat[,2]<-jitter(tl$TL,0.1)              # y-axis value based on trophic level
  lMat[,1]<-runif(vcount(g))               # randomly assign along x-axis

  colTL <-as.numeric(cut(tl$TL,11))   # Divide trophic levels in 11 segments
  colnet <- brewer.pal(11,"RdYlGn")   # Assign colors to trophic levels
  V(g)$color <- colnet[12-colTL]      # Add colors to the igraph object


  if(tk){
    tkid <- tkplot(g, edge.width=.3,edge.arrow.size=.4,
       vertex.label.color="white",
       edge.color="grey50",
       edge.curved=0.3, layout=lMat)
    return( tkplot.getcoords(tkid))

  } else {
    plot(g, edge.width=.3,edge.arrow.size=.4,
         vertex.label.color="white",
         edge.color="grey50",
         edge.curved=0.3, layout=lMat)

  }


}



#' Calculate topological indices for ecological networks
#'
#' @param ig vector of igraph objects
#'
#' @return a data.frame with the following fields:
#'
#'         \item{Size:}{Number of species}
#'         \item{Top:}{Number of top predator species}
#'         \item{Basal:}{Number of basal especies}
#'         \item{Omnivory:}{ Number of omnivory species, is the number }
#'         \item{Links:}{ number of interactions}
#'         \item{LD:}{ linkage density}
#'         \item{Connectance:}{ Connectance}
#'         \item{PathLength:}{ average path length}
#'         \item{Clustering:}{ clustering coeficient}
#'         \item{Cannib:}{ number of cannibalistic species}
#'         \item{TLmean:}{ mean trophic level}
#'         \item{Cannib:}{ maximum trophic level}
#'
#'
#' @export
#' @importFrom NetIndices TrophInd
#' @import igraph
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

topologicalIndices <- function(ig){

  if(class(ig)=="igraph") {
    ig <- list(ig)
  } else if(class(ig[[1]])!="igraph") {
    stop("parameter ig must be an igraph object")
  }

  df <- lapply(ig, function(g){

    deg <- degree(simplify(g), mode="out") # calculate the out-degree: the number of predators

    V(g)$outdegree <-  deg

    nTop <- length(V(g)[outdegree==0]) # Top predators do not have predators

    deg <- degree(g, mode="in") # calculate the in-degree: the number of preys

    V(g)$indegree <-  deg

    nBasal <- length(V(g)[indegree==0]) # Basal species do not have preys

    vcount(g)-nTop-nBasal

    size <- vcount(g)

    links <- ecount(g)

    linkDen <- links/size          # Linkage density

    conn <- links/size^2           # Connectance

    pathLength <- average.path.length(g)   # Average path length

    clusCoef <- transitivity(g, type = "global")

    cannib <- sum(which_loop(g))

    tl <- TrophInd(get.adjacency(g,sparse=F))  # Calculate the trophic level

    omn <- sum(round(tl$OI,5)>0)

    data.frame(Size=size,Top=nTop,Basal=nBasal,Omnivory=omn,Links=links, LD=linkDen,Connectance=conn,PathLength=pathLength,
      Clustering=clusCoef, Cannib=cannib, TLmean=mean(tl$TL),TLmax=max(tl$TL))
  })

  do.call(rbind,df)
}


#' Paralell version of topologicalIndices
#'
#' @param ig vector of igraph objects
#'
#' @return the same as topologicalIndices
#' @export
#'
#' @examples
#'
#' parTopologicalIndices(netData)
parTopologicalIndices <- function(ig){

  cn <-detectCores()
  cl <- makeCluster(cn)
  registerDoParallel(cl)

  df <- foreach(i=seq_along(ig), .combine='rbind',.inorder=FALSE,.packages='igraph') %dopar% {

    topologicalIndices(g)

  }
  stopCluster(cl)

  return(df)
}


#' Calculate the incoherence index of a food web
#'
#' The incoherence index is based in how the species fit in discrete trophic levels
#' when Q is closer to 0 more coherent and stable is a food web.
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
#'
#' @return a data.frame with the following fields
#'
#'       \item{Q}{incoherence (0=coherent)}
#'      \item{rQ}{ratio of Q with expected Q under null expectation of a random network given N=nodes L=links B=basal nodes Lb=basal links}
#'      \item{mTI}{mean trophic level}
#'      \item{rTI}{ratio of mTI with expected TI under the same null model expectation}
#
#'
#' @export
#'
#' @examples
#'
#' calcIncoherence(netData[[1]])
#'
#'
calcIncoherence <- function(g,ti=NULL) {
  require(igraph)
  require(NetIndices)
  if(is.null(ti) )
    ti<-TrophInd(get.adjacency(g,sparse=FALSE))
  v <- ti$TL
  z <- outer(v,v,'-');
  A <- get.adjacency(g,sparse = FALSE)
  xx <- A>0
  x <- (A*t(z))[xx]
  #meanQ <- sum(x)/ecount(g)
  #sdQ <- sqrt(sum((x-1)^2)/vcount(g) )
  Q <- sqrt(sum(x*x-1)/ecount(g) )

  basal <- which(round(v,8)==1)
  bedges <- sum(degree(g,basal,mode='out'))
  mTI <- mean(v)
  mK <- mean(degree(g,mode='out'))
  eTI <- 1+(1-length(basal)/vcount(g))*ecount(g)/bedges
  eQ <- sqrt(ecount(g)/bedges-1)
  data.frame(Q=Q,rQ=Q/eQ,mTI=mTI,rTI=mTI/eTI)
}
