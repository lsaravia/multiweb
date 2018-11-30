
#' Plot ecological network organized by trophic level, with node size determined by the node degree
#'
#' @param ig igraph object
#' @param vertexLabel logical plot vertex labels
#' @param vertexSizeFactor numeric factor to determine the size of the label with degree
#' @param tk TRUE generate an interactive plot using tkplot
#' @param modules FALSE plot modules in the x axis, obtained with the cluster spinglass algorithm
#' @param lMat Matrix of postions for the nodes
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
plotTrophLevel <- function(g,vertexLabel=FALSE,vertexSizeFactor=5,tk=FALSE,modules=FALSE,lMat=NULL, ...){

  deg <- degree(g, mode="all") # calculate the degree: the number of edges
  # or interactions

  V(g)$size <- log10(deg)*vertexSizeFactor+vertexSizeFactor    # add node degrees to igraph object

  V(g)$frame.color <- "white"    # Specify plot options directly on the object

  V(g)$color <- "orange"         #
  E(g)$color <- "gray50"

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
  if(is.null(lMat)){

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
          if(length(m$membership)==0)
            m$membership <- 1
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
  }

  colTL <-as.numeric(cut(tl$TL,11))   # Divide trophic levels in 11 segments
  colnet <- brewer.pal(11,"RdYlGn")   # Assign colors to trophic levels
  V(g)$color <- colnet[12-colTL]      # Add colors to the igraph object


  if(tk){
    tkid <- tkplot(g, edge.width=.3,edge.arrow.size=.4,
       vertex.label.color="white",
       edge.curved=0.3, layout=lMat)
    return( tkplot.getcoords(tkid))

  } else {
    plot(g, edge.width=.3,edge.arrow.size=.4,
         vertex.label.color="white",
         edge.curved=0.3, layout=lMat,...)

  }

}


