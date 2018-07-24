# Function definitions



#' Read ecological networks in CSV  or tab separated file format as edge list or adyacency matrix
#'
#' @param fileName vector of fileNames with the networks
#' @param filePath path of the files NULL by default
#'
#' @return an igraph object if there is only one file or a list of igraph objects named after the list without extension
#' @export
#'
#' @examples
#'
#' # Reads a network in edge list (interaction list) format, with predators as the first column
#' #
#' fileName <- system.file("extdata", "WeddellSea_FW.csv", package = "EcoNetwork"))
#' g <- readEcoNetwork(fileName)
#'
#' # Reads a network in adyacency matrix format, with predators as columns
#' #
#' fileName <- system.file("extdata", "BarentsBoreal_FW.csv", package = "EcoNetwork"))
#' g <- readEcoNetwork(fileName)
#'
#' # Read a vector of files
#' #
#'\dontrun{
#' dn <- list.files("inst/extdata",pattern = "^.*\\.csv$")
#' netData <- readEcoNetwork(dn,"inst/extdata")
#'}

readEcoNetwork <- function(fileName,filePath=NULL){

  fn <-  if(!is.null(filePath)) paste0(filePath,"/",fileName) else fileName

  g <- lapply(fn, function(fname){

    fe <- tools::file_ext(fname)
    if(fe=="csv")
      web <- read.csv(fname,  header = T,check.names = F,stringsAsFactors = FALSE)
    else
      web <- read.delim(fname,  header = T,check.names = F,stringsAsFactors = FALSE)

    if( ncol(web)==2 ){
      web <- web[,c(2,1)]

      g <- igraph::graph_from_data_frame(web)

    } else {
      if( (ncol(web)-1) == nrow(web)  ) {                   # The adjacency matrix must be square
        g <- igraph::graph_from_adjacency_matrix(as.matrix(web[,2:ncol(web)]))

      } else {
        g <- NULL
        warning("Invalid file format: ",fileName)
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
#'
#' @return a plot
#' @export
#'
#'
#' @importFrom NetIndices TrophInd
#' @importFrom igraph     V degree
#' @importFrom RColorBrewer brewer.pal
#' @examples
#'
#' plotEcoNetworkTrophLevel(netData[[1]])
plotEcoNetworkTrophLevel <- function(g,vertexLabel=FALSE,vertexSizeFactor=5){

  deg <- degree(g, mode="all") # calculate the degree: the number of edges
  # or interactions

  V(g)$size <- log10(deg)*vertexSizeFactor+vertexSizeFactor    # add node degrees to igraph object

  V(g)$frame.color <- "white"    # Specify plot options directly on the object

  V(g)$color <- "orange"         #

  if(!vertexLabel)
    V(g)$label <- NA

  tl <- NetIndices::TrophInd(get.adjacency(g,sparse=F))  # Calculate the trophic level

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


  plot(g, edge.width=.3,edge.arrow.size=.4,
       #vertex.label=NA,
       vertex.label.color="white",
       edge.color="grey50",
       edge.curved=0.3, layout=lMat)


}



#' Calculate topological indices for ecological networks
#'
#' @param ig vector of igraph objects
#'
#' @return a data.frame with the following fields:
#'         Size: Number of species
#'         Top:  Number of top predator species
#'         Basal: Number of basal especies
#'         Links: number of interactions
#'         LD: linkage density
#'         Connectance: Connectance
#'         PathLength: average path length
#'         Clustering: clustering coeficient
#'         Cannib: number of cannibalistic species
#'
#' @export
#'
#' @examples
topologicalIndicesEcoNetwork <- function(ig){

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

    data.frame(Size=size,Top=nTop,Basal=nBasal,Links=links, LD=linkDen,Connectance=conn,PathLength=pathLength,Clustering=clusCoef, Cannib=cannib)
  })

  do.call(rbind,df)
}


parTopologicalIndicesEcoNetwork <- function(ig){

  cn <-detectCores()
  cl <- makeCluster(cn)
  registerDoParallel(cl)

  df <- foreach(i=seq_along(ig), .combine='rbind',.inorder=FALSE,.packages='igraph') %dopar% {

    g <- ig[[i]]

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

    data.frame(Size=size,Top=nTop,Basal=nBasal,Links=links, LD=linkDen,Connectance=conn,PathLength=pathLength,Clustering=clusCoef, Cannib=cannib)
  }
  stopCluster(cl)

  return(df)
}
