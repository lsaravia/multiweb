# Function definitions



#' Read ecological networks in CSV or tab separated file format as edge list or adyacency matrix
#'
#' @param fileName vector of fileNames with the networks
#' @param filePath path of the files NULL by default
#' @param fhead TRUE if the files have header fields, FALSE otherwise.
#' @param skipColum integer, number of columns that are skiped 1 by default
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

readNetwork <- function(fileName,filePath=NULL,fhead=TRUE,skipColumn=1){

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
      if( (ncol(web)-skipColumn) == nrow(web)  ) {                   # The adjacency matrix must be square
        skipColumn <- skipColumn+1
        g <- igraph::graph_from_adjacency_matrix(as.matrix(web[,skipColumn:ncol(web)]),mode="directed")

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

#' Read ecological multiplex networks using different files for each layer as a 'node-aligned' multilayer network
#'
#' This functions uses [readNetwork()] to read files that represent network layers that can be
#'   different interaction types, represented by the type attribute of the igraph object.
#'
#' @param fileName vector of fileNames with the layers of the networks
#' @param types vector of types that represent the layers
#' @param filePath path of the files NULL by default
#' @param fhead TRUE if the files have header fields, FALSE otherwise.
#' @param skipColum integer, number of columns that are skiped 1 by default
#' @param format string, "layers" is the default were diferent layers are coded as different files
#'   there must be the same number of files as the length of the types vector as each type represent
#'   a layer. "GLV" Represent multiple intarction types as pairs of entries in a matrix so competition
#'   is represented as a[i,j]= -1, a[j,i]=-1, predation a[i,j]=1,a[j,i]=-1 where species j is
#'   the predator. Mutualism is a[i,j]=a[j,i]=1. Any negative or positive number works because intensity of interactions are not registered.
#'
#' @return an mgraph object inherited from igraph object  with the type attribute set to each kind of layer or interaction
#' @export
#'
#' @importFrom igraph     E graph_from_data_frame graph_from_adjacency_matrix V
#' @seealso [readNetwork()]
#' @examples
#'
#' # Read a vector of files
#' #
#'\dontrun{
#'
#' fpath <- system.file("extdata", package = "EcoNetwork")
#' dn <- list.files(fpath,pattern = "^Kefi2015.*\\.txt$")
#' netData <- readMultiplex(dn,c("Competitive","Mutualistic","Trophic"),fpath,skipColum=2)
#'}

readMultiplex <- function(fileName,types=c("Competitive","Mutualistic","Trophic"),filePath=NULL,
                          fhead=TRUE,skipColumn=1,format="layers"){

  if( format=="layers"){
    g <- readNetwork(fileName=fileName,filePath=filePath, fhead=fhead,skipColumn = skipColumn)
    stopifnot(length(g)==length(types))

    for(t in seq_along(types)){
      E(g[[t]])$type <- types[t]
    }
    gg <- lapply(g, igraph::as_data_frame)
    gt <- do.call(rbind,gg)
    gt <- graph_from_data_frame(gt)

  } else if(format=="GLV"){

    stopifnot(length(types)==3)

    fn <-  if(!is.null(filePath)) paste0(filePath,"/",fileName) else fileName

    web <- read.delim(fn,stringsAsFactors = FALSE)
    if( (ncol(web)-skipColumn) == nrow(web)  ) {                 # The adjacency matrix must be square
      skipColumn <- skipColumn+1
      web <- web[,skipColumn:ncol(web)]
      pred <- matrix(0,nrow = nrow(web),ncol = nrow(web))
      comp <- matrix(0,nrow = nrow(web),ncol = nrow(web))
      mut  <- matrix(0,nrow = nrow(web),ncol = nrow(web))

      for( i in 1:nrow(web))
        for(j in 1:nrow(web)){
          if(web[i,j]!=0){
            if(web[i,j]>0 && web[j,i]<0){
              pred[i,j]<-1
            } else if(web[i,j]<0 && web[j,i]==0) {
              comp[i,j] <- 1
            }else if(web[i,j]<0 && web[j,i]<0) {
              comp[i,j] <- 1
              comp[j,i] <- 1
            } else if(web[i,j]>0 && web[j,i]==0) {
              mut[i,j] <- 1
            } else if(web[i,j]>0 && web[j,i]>0){
              mut[i,j] <- 1
              mut[j,i] <- 1
            } else if(web[i,j]<0 && web[j,i]>0){
              pred[j,i]<-1
            } else {
              warning("Unclasified interaction type ", web[i,j], web[j,i])
            }

          }

        }

      gg <- list(comp,mut,pred)
      # convert to igraph
      gdf <- lapply(seq_along(gg), function(t){
        g <- graph_from_adjacency_matrix(gg[[t]])
        E(g)$type <- types[t]
        df <- igraph::as_data_frame(g)
      })
      gt <- do.call(rbind,gdf)
      gt <- graph_from_data_frame(gt)
      V(gt)$name <- names(web)

    } else {
      gt <- NULL
      warning("Interaction matrix not square")
    }

  }
  class(gt) <- c(class(gt),"mgraph")

  return(gt)
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
#' @return if tk==TRUE returns a layout matrix, else returns a plot
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

  if(inherits(g, "mgraph") && ("Trophic" %in% unique(unlist(edge.attributes(g)))) ){
    tt <- subgraph.edges(g,E(g)[E(g)$type=="Trophic"])
    tl <- TrophInd(get.adjacency(tt,sparse=F))
  } else {
    tl <- TrophInd(get.adjacency(g,sparse=F))  # Calculate the trophic level
  }
  # Layout matrix to specify the position of each vertix
  # Rows equal to the number of vertices (species)
  lMat <-matrix(
    nrow=vcount(g),
    ncol=2
  )

  lMat[,2]<-jitter(tl$TL,0.1)              # y-axis value based on trophic level

  if(modules) {
    if(count_components(g)>1){
      if(!is.named(g)) V(g)$name <- (1:vcount(g))
      dg <- components(g)
      V(g)$membership = 0
      for(comp in unique(dg$membership)) {
        g1 <- induced_subgraph(g, which(dg$membership == comp))
        m<-cluster_spinglass(g1)
        V(g)[V(g1)$name]$membership <-  m$membership + max(V(g)$membership)
      }
      m$membership <- V(g)$membership

    } else {
      m<-cluster_spinglass(g)
    }

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


