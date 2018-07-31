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
#' @param modules FALSE plot modules in the x axis, obtained with the cluster spinglass algorithm
#' @param title title for the plot
#' @param ... Addittional parameters to the plot function
#'
#' @return a plot
#' @export
#'
#'
#' @importFrom NetIndices TrophInd
#' @importFrom igraph     V degree get.adjacency cluster_spinglass
#' @importFrom RColorBrewer brewer.pal
#' @examples
#'
#' plotTrophLevel(netData[[1]])
plotTrophLevel <- function(g,vertexLabel=FALSE,vertexSizeFactor=5,tk=FALSE,modules=FALSE, ...){

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

  if(modules) {
    m<-cluster_spinglass(g)
    lMat[,1]<-jitter(m$membership,1) # randomly assign along x-axis
  } else {
    lMat[,1]<-runif(vcount(g))               # randomly assign along x-axis
  }

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
         edge.curved=0.3, layout=lMat,...)

  }

}


