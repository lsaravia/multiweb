
#' Curbe ball algorithm to generate random networks with given an igraph network
#'
#' This extract the adjacency matrix from an igraph object and generates a randomized version
#' of the network with the same row and column totals, the results have the same in and out degree sequence
#' than the original network. The diagonals are not avoided so it can generate self-links or cannibalism in
#' the context of food-webs. In the case that the network has multiple components (or disconnected networks)
#' the algorithm simulates around the same number of components, if the original network has 1 component
#' the algorithm enforces that the results have all 1 component.
#'
#' Based on:
#'
#' @references Strona, G. et al. 2014. A fast and unbiased procedure to randomize ecological binary matrices with
#' fixed row and column totals. -Nat. Comm. 5: 4114. doi: 10.1038/ncomms5114
#'
#' @param g igraph object to extract adjacency matrix
#' @param nsim number of generated random networks
#'
#' @return a list of randomized igraph objects
#' @export
#'
#' @examples
#'
#' curveBall(netData[[1]])
curveBall<-function(g,nsim=1000){
  stopifnot(class(g)=="igraph")
  m <- get.adjacency(g,sparse=FALSE)

  if(components(g)$no>1) {

    nets<- lapply(1:nsim, function (x) {
      # repeat {
      RC <- dim(m)
      R  <- RC[1]
      C  <- RC[2]
      hp <- list()
      for (row in 1:dim(m)[1]) {hp[[row]] <- (which(m[row,]==1))}
      l_hp <- length(hp)
      for (rep in 1:5*l_hp){
        AB <- sample(1:l_hp,2)
        a  <- hp[[AB[1]]]
        b  <- hp[[AB[2]]]
        ab <- intersect(a,b)
        l_ab <- length(ab)
        l_a <- length(a)
        l_b <- length(b)
        if ((l_ab %in% c(l_a,l_b))==FALSE){
          tot <- setdiff(c(a,b),ab)
          l_tot <- length(tot)
          tot <- sample(tot, l_tot, replace = FALSE, prob = NULL)
          L <- l_a-l_ab
          hp[[AB[1]]]  <-  c(ab,tot[1:L])
          hp[[AB[2]]]  <-  c(ab,tot[(L+1):l_tot])}

      }
      rm <- matrix(0,R,C)
      for (row in 1:R){rm[row,hp[[row]]] <- 1}

      g <- graph_from_adjacency_matrix(rm,mode="directed")
      return(g)
    })
  } else {   #if the the networks has only one component enforce that in the simulations

    nets<- lapply(1:nsim, function (x) {
      repeat {
      RC <- dim(m)
      R  <- RC[1]
      C  <- RC[2]
      hp <- list()
      for (row in 1:dim(m)[1]) {hp[[row]] <- (which(m[row,]==1))}
      l_hp <- length(hp)
      for (rep in 1:5*l_hp){
        AB <- sample(1:l_hp,2)
        a  <- hp[[AB[1]]]
        b  <- hp[[AB[2]]]
        ab <- intersect(a,b)
        l_ab <- length(ab)
        l_a <- length(a)
        l_b <- length(b)
        if ((l_ab %in% c(l_a,l_b))==FALSE){
          tot <- setdiff(c(a,b),ab)
          l_tot <- length(tot)
          tot <- sample(tot, l_tot, replace = FALSE, prob = NULL)
          L <- l_a-l_ab
          hp[[AB[1]]]  <-  c(ab,tot[1:L])
          hp[[AB[2]]]  <-  c(ab,tot[(L+1):l_tot])}

      }
      rm <- matrix(0,R,C)
      for (row in 1:R){rm[row,hp[[row]]] <- 1}

      g <- graph_from_adjacency_matrix(rm,mode="directed")
      if(components(g)$no==1)
         break
      }
      return(g)
    })
  }
}

#' Generate directed Erdos-Renyi random networks with at least 1 basal node and only one component
#'
#' This uses the igraph's function sample_gnm to generate nsim random networks with the same number of nodes
#' and links than the parameter ig and two restrictions:
#' 1) at least one basal species/node, that is a species that has no prey, 2) 1 connected component so there is no
#' disconnected species or sub-community.
#'
#' @param ig igraph object with parameters to use in the random network simulations: number of species/nodes
#'           and number of links/edges
#' @param nsim number of simulations
#'
#' @return a list with igraph objects
#' @export
#'
#' @examples
#'
#' generateERbasal(netData[[1]])
generateERbasal <- function(ig,nsim=1000){
  if(!is_igraph(ig))
    stop("Parameter ig must be an igraph object")

  size <- vcount(ig)

  links <- ecount(ig)

  er <- lapply(1:nsim, function (x) {
    e <- sample_gnm(size, links, directed = TRUE)
    basal <- length(V(e)[degree(e,mode="in")==0])
    while(components(e)$no>1 | basal==0){
      e <- sample_gnm(size, links,directed = TRUE)
      basal <- length(V(e)[degree(e,mode="in")==0])
    }
    return(e) })
}


