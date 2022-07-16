# Function definitions



#' Read ecological networks in CSV or tab separated file format as edge list or adyacency matrix
#'
#' If the network is in edge list format and there is a third column is treated as an edge attribute with the name of the column
#'
#' @param fileName vector of fileNames with the networks
#' @param filePath path of the files NULL by default
#' @param fhead TRUE if the files have header fields, FALSE otherwise.
#' @param skipColum integer, number of columns that are skiped 1 by default
#' @param edgeListFormat integer, for the edge list format, if 1 the first column is the "in" (the predator),
#'                       if 2 the first column is the "out" link (the prey).
#'
#' @return an igraph object if there is only one file or a list of igraph objects named after the list without extension
#' @export
#'
#' @examples
#'
#' # Reads a network in edge list (interaction list) format, with predators as the first column by default
#' #
#' fileName <- system.file("extdata", "WeddellSea_FW.csv", package = "multiweb")
#' g <- readNetwork(fileName)
#'
#' # Reads a network in adyacency matrix format, with predators as columns
#' #
#' fileName <- system.file("extdata", "BarentsBoreal_FW.csv", package = "multiweb")
#' g <- readNetwork(fileName)
#'
#' # Read a vector of files
#' #
#'\dontrun{
#' dn <- list.files("inst/extdata",pattern = "^.*\\.csv$")
#' netData <- readNetwork(dn,"inst/extdata")
#'}

readNetwork <- function(fileName,filePath=NULL,fhead=TRUE,skipColumn=1,edgeListFormat=1){

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

    if( ncol(web)<=3 ){
      if( edgeListFormat==1 )
        if( ncol(web)==2)
          web <- web[,c(2,1)]
        else
          web <- web[,c(2,1,3)]

      g <- igraph::graph_from_data_frame(web)                        # the 3d field generate an attribute

    } else {
      if( (ncol(web)-skipColumn) == nrow(web)  ) {                   # The adjacency matrix must be square
        skipColumn <- skipColumn+1
        g <- igraph::graph_from_adjacency_matrix(as.matrix(web[,skipColumn:ncol(web)]),mode="directed")
        if(skipColumn>1){
          #
          # Assume first column before data is species names
          #
          nmes <- web[,skipColumn-1]
          igraph::V(g)$name <- nmes
        } else {
          # if all col names starts with V strip the V
          #
          if( all(grepl("^V",names(web))) ) {
            names(web) <- sub("^V(*.)","\\1",names(web))
            igraph::V(g)$name <- names(web)
            }
        }
        g
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
#' @return an mgraph object that is a list of igraph objects with the type attribute set to each kind of layer or interaction
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
#' fpath <- system.file("extdata", package = "multiweb")
#' dn <- list.files(fpath,pattern = "^Kefi2015.*\\.txt$")
#' netData <- readMultiplex(dn,c("Competitive","Mutualistic","Trophic"),fpath,skipColum=2)
#'}

readMultiplex <- function(fileName,types=c("Competitive","Mutualistic","Trophic"),filePath=NULL,
                          fhead=TRUE,skipColumn=1,format="layers"){

  if( format=="layers"){
    gt <- readNetwork(fileName=fileName,filePath=filePath, fhead=fhead,skipColumn = skipColumn)
    stopifnot(length(gt)==length(types))

    for(t in seq_along(types)){
      E(gt[[t]])$type <- types[t]
    }
    names(gt) <- types

#    gg <- lapply(g, igraph::as_data_frame)
#    gt <- do.call(rbind,gg)
#    gt <- graph_from_data_frame(gt)

  } else if(format=="GLV"){

    stopifnot(length(types)==3)

    fn <-  if(!is.null(filePath)) paste0(filePath,"/",fileName) else fileName

    web <- read.delim(fn,stringsAsFactors = FALSE)
    if( (ncol(web)-skipColumn) == nrow(web)  ) {                 # The adjacency matrix must be square

      nmes <- web[,skipColumn]
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
      gt <- lapply(seq_along(gg), function(t){
        g <- graph_from_adjacency_matrix(gg[[t]])
        E(g)$type <- types[t]
        V(g)$name <- nmes
        g
        #df <- igraph::as_data_frame(g)
      })
      #gt <- do.call(rbind,gdf)
      #gt <- graph_from_data_frame(gt)
      #V(gt)$name <- nmes
      names(gt) <- types

    } else {
      gt <- NULL
      warning("Interaction matrix not square")
    }

  }
  class(gt) <- c("mgraph")

  return(gt)
}



#' From multiple interaction object 'mgraph'to GLV adjacency matrix
#'
#' This functions takes a 'mgraph' object and convert it to a Generalized Lotka-Volterra adjacency matrix
#' Position is important in fact the order is  Competitive/Negative,Mutualistic/Positive,Trophic/Antagonistic
#' If 'istrength' is TRUE the attribute weight is assumed as the strength of the interacion
#'
#'
#' @param mg multiple interaction object, class 'mgraph'
#' @param types vector of types that represent the layers
#' @param istrength edge weights
#'
#' @return an matrix with the signs of the interactions
#' @export
#'
#' @importFrom igraph     E graph_from_data_frame graph_from_adjacency_matrix V
#' @seealso [readMultiplex()]
#' @examples
#'
#' # Read a vector of files
#' #
#'\dontrun{
#'
#' fpath <- system.file("extdata", package = "multiweb")
#' dn <- list.files(fpath,pattern = "^Kefi2015.*\\.txt$")
#' netData <- readMultiplex(dn,c("Competitive","Mutualistic","Trophic"),fpath,skipColum=2)
#' toGLVadjMat(netData)
#'}

toGLVadjMat <- function(mg,types=c("Competitive","Mutualistic","Trophic"),istrength=FALSE){

  if( class(mg)!='mgraph')
    stop("parameter mg must be an mgraph object")
  if( length(mg)!=3 )
    stop("parameter mg must have 3 components")

  stopifnot(length(types)==3)

  if(is.null(edge_attr(mg[[types[3]]],"weight")))
    pred <- as_adj(mg[[types[3]]],sparse = FALSE)
  else
    pred <- as_adj(mg[[types[3]]],sparse = FALSE,attr="weight")

  if(is.null(edge_attr(mg[[types[1]]],"weight")))
    comp <- as_adj(mg[[types[1]]],sparse = FALSE)
  else
    comp <- as_adj(mg[[types[1]]],sparse = FALSE,attr="weight")

  if(is.null(edge_attr(mg[[types[2]]],"weight")))
    mut <- as_adj(mg[[types[2]]],sparse = FALSE)
  else
    mut <- as_adj(mg[[types[2]]],sparse = FALSE,attr="weight")

  web  <- matrix(0,nrow = nrow(pred),ncol = nrow(pred))

  if(istrength)
  {
    for( i in 1:nrow(web)){
      for(j in 1:nrow(web)){
        if(pred[i,j]){
            web[i,j] <- -pred[i,j]
            web[j,i] <- pred[i,j]
        } else if(comp[i,j]) {

            web[i,j] <- -comp[i,j]

        } else if(mut[i,j]) {
            web[i,j] <- mut[i,j]
        }

      }

    }

  } else {

    for( i in 1:nrow(web)){
      for(j in 1:nrow(web)){
        if(pred[i,j]){
            web[i,j] <- -1
            web[j,i] <- 1
        } else if(comp[i,j]) {

            web[i,j] <- -1

        } else if(mut[i,j]) {
            web[i,j] <- 1
        }

      }

    }
  }

  return(web)
}



#' From multiple igraph objects to 'mgraph'
#'
#' This functions takes a list of igraph objects and convert it to a 'mgraph' object.
#'
#'
#' @param mg list of 'igraph' objects
#' @param types vector of types that represent the layers with the same length as mg
#'
#' @return class 'mgraph' object
#' @export
#'
#' @seealso [readMultiplex()]
#' @examples
#'
#' # Read a vector of files
#' #
#'\dontrun{
#'
#' fileName <- c(system.file("extdata",  package = "multiweb"))
#' dn <- list.files("inst/extdata",pattern = "^Kefi2015.*\\.txt$")
#' g <- readNetwork(dn,"inst/extdata", skipColumn = 2)
#' fromIgraphToMgraph(g,c("Negative","Positive","Antagonistic"))
#'
#'}

fromIgraphToMgraph <- function(g,types){

  if( any(sapply(g,class)!='igraph'))
    stop("parameter mg must be a list of igraph objects")

  if(length(g)!= length(types))
     stop("Length of mg must be equal to length of types")

  names(g) <- types
  class(g) <- 'mgraph'

  return(g)
}

#' From Generalized Lotka Volterra adjacency matrix to igraph topological object
#'
#' It counts predator-prey/Antagonistic interactions like 1 edge,
#' competition and mutualisms are counted as is: two edges or one edge
#'
#' @param glvAdj numeric matrix Generalized Lotka-Volterra adjacency matrix
#' @param spc numeric vector of present Species
#'
#' @return an igraph object
#' @export
#' @importFrom igraph     graph_from_adjacency_matrix V
#'
#' @examples
#' # Build a matrix
#'
#' m <- matrix(0,nrow=5,ncol=5)
#' m[1,2] <- m[1,3] <- m[3,4]<- .2
#' m[2,1] <- m[3,1] <- m[4,3] <- -2
#' m[5,4] <- m[4,5] <- 0.1            # Mutualistic
#' m[1,1] <- -0.01                    # Cannibalistic
#'
#' g <- fromGLVadjToIgraph(m,c(1,1,1,1,0))
#'
#' g <- fromGLVadjToIgraph(m,c(0,1,1,1,1))
#'
fromGLVadjToIgraph<- function(glvAdj,spc){

  stopifnot(nrow(glvAdj)==ncol(glvAdj))
  if(length(spc)!=nrow(glvAdj)) stop("length of species vector has to be equal to glvAdj dimensions")

  d <- spc!=0
  A <- glvAdj[d,d]
  for(i in seq_len(nrow(A)))
    for(j in seq_len(nrow(A))){
      if(A[i,j]<0 && A[j,i]>0){
        A[i,j] <- 1
        A[j,i] <- 0
      }
    }
  # diag(A) <- 0
  A[A>0] <- 1
  A[A<0] <- 1
  g <- graph_from_adjacency_matrix(A,mode="directed")
  V(g)$name <- which(spc>0)                              # Maintains the "names" of the original matrix
  return(g)
}
