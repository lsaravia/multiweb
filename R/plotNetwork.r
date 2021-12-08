
#' Plot ecological network organized by trophic level, with node size determined by the node degree, and modules.
#'
#' @param ig igraph object
#' @param vertexLabel logical plot vertex labels
#' @param vertexSizeFactor numeric factor to determine the size of the label with degree
#' @param vertexSizeMin    numeric determine the minimum size of the label
#' @param tk TRUE generate an interactive plot using tkplot and returns a matrix with coordinates from [igraph::tkplot.getcoords()]
#' @param lMat Matrix of postions for the nodes
#' @param modules if TRUE plot modules in the x axis, obtained with the cluster spinglass algorithm
#' @param weights The weights of the edges for the [igraph::cluster_spinglass()] community detection, either a numeric vector, NULL, or NA.
#'                if NULL and the network has a 'weight' attribute then that will be used,
#'                if NA then the 'weight' attributte is not considered.
#' @param community_obj Insteado of calculating modules/communities with cluster spinglass take a community object
#' @param bpal if NULL it uses the "RdYlGn" RColorBrewer palette, else must be a vector of colors of length 11.
#' @param maxTL maximum trophic level to draw y-axis
#'
#' @param ... Addittional parameters to the plot function
#'
#' @return returns a plot and if tk==TRUE returns a layout matrix
#' @export
#'
#' @aliases plotTrophLevel
#'
#'
#' @importFrom NetIndices TrophInd
#' @importFrom igraph     V degree get.adjacency cluster_spinglass vcount
#' @importFrom RColorBrewer brewer.pal
#' @importFrom dplyr      "%>%" mutate dense_rank
#' @examples
#'
#' plot_troph_level(netData[[1]])
plot_troph_level <- function(g,vertexLabel=FALSE,vertexSizeFactor=5,vertexSizeMin=5,tk=FALSE,modules=FALSE,lMat=NULL,weights=NA,community_obj=NULL, bpal= NULL,
                             maxTL=NULL,...){

  deg <- degree(g, mode="all") # calculate the degree: the number of edges
  # or interactions

  V(g)$size <- log10(deg)*vertexSizeFactor+vertexSizeMin    # add node degrees to igraph object

  V(g)$frame.color <- "white"    # Specify plot options directly on the object

  V(g)$color <- "orange"         #
  E(g)$color <- "gray50"
  #E(g)$width <- .3
  # if( !is.null(edge.width )) {
  #   if( is.null(weights) )
  #       E(g)$width <- E(g)$width* edge.width
  #   else if( is.vector(weights, mode="numeric"))
  #       E(g)$width <- weights* edge.width
  # }

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
      if(!is.null(community_obj)) {
        m <- community_obj
      } else {
        if(count_components(g)>1){
          if(!is.named(g)) V(g)$name <- (1:vcount(g))
          dg <- components(g)
          V(g)$membership = 0
          for(comp in unique(dg$membership)) {
            g1 <- induced_subgraph(g, which(dg$membership == comp))
            m<-cluster_spinglass(g1,weights=weights)
            if(length(m$membership)==0)
              m$membership <- 1
            V(g)[V(g1)$name]$membership <-  m$membership + max(V(g)$membership)
          }
          m$membership <- V(g)$membership

        } else {
          m<-cluster_spinglass(g,weights=weights)
        }
      }
      # Order groups in ascending trofic level
      #
      df <- data.frame(tl=tl$TL,m=m$membership)
      df <- df %>% mutate(m = dense_rank(ave(tl, m, FUN = max)))
      lMat[,1]<-jitter(df$m,1) # randomly assign along x-axis

    } else {
      lMat[,1]<-runif(vcount(g))               # randomly assign along x-axis
    }
  }

  colTL <-as.numeric(cut(tl$TL,11))   # Divide trophic levels in 11 segments
  if( !is.null(bpal)) {
    colnet <- colorRampPalette(bpal)(11)
  } else {
    colnet <- brewer.pal(11,"RdYlGn")   # Assign colors to trophic levels
  }
  V(g)$color <- colnet[12-colTL]      # Add colors to the igraph object


  if(tk){
    tkid <- tkplot(g, edge.arrow.size=.4,
       edge.curved=0.3, layout=lMat,...)
    return( tkplot.getcoords(tkid))

  } else {
    plot(g, edge.arrow.size=.4,
         edge.curved=0.3, layout=lMat,...)
    maxnew <- max(tl$TL)
    minnew <- min(tl$TL)
    maxold <- 1
    minold <- -1
    t2 <- function(x) (maxold-minold)/(maxnew -minnew)*(x - maxnew)+maxold
    tlseq <- seq(1,ifelse(is.null(maxTL),maxnew+1,maxTL),by=1)
    axis(side=2,at=t2(tlseq),labels=tlseq,  las=1, col = NA, col.ticks = 1)


  }

}


#' @export
#'
plotTrophLevel <- function(g,vertexLabel=FALSE,vertexSizeFactor=5,vertexSizeMin=5,tk=FALSE,modules=FALSE,lMat=NULL,weights=NA,community_obj=NULL, bpal= NULL,
                           maxTL=NULL,...)
  {
  plot_troph_level(g,vertexLabel,vertexSizeFactor,vertexSizeMin,tk,modules,lMat,weights, community_obj,bpal,maxTL, ...)
}


