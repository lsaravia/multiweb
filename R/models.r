
#' Curbe ball algorithm to generate random networks with given an igraph network
#'
#' This extract the adjacency matrix from an igraph object and generates a randomized version
#' of the network with the same row and column totals, the results have the same in and out degree sequence
#' than the original network. The diagonals are not avoided so it can generate self-links or cannibalism in
#' the context of food-webs. In the case that the network has multiple components (or disconnected networks)
#' the algorithm simulates around the same number of components, if the original network has 1 component
#' the algorithm enforces that the results have all 1 component. If the edge attribute weight is present
#' and istrength=TRUE then the weigth is additionally randomized keeping the column sum equal to the original matrix.
#'
#' Based on:
#'
#' @references Strona, G. et al. 2014. A fast and unbiased procedure to randomize ecological binary matrices with
#' fixed row and column totals. -Nat. Comm. 5: 4114. doi: 10.1038/ncomms5114
#'
#' @aliases curveBall
#' @param g igraph object to extract adjacency matrix
#' @param nsim number of generated random networks
#' @param istrength if TRUE the edge attribute weight is taken as interaction strength and randomized keeping the column sum equal to the original matrix.
#'
#' @return a list of randomized igraph objects
#' @export
#'
#' @examples
#'
#' curve_ball(netData[[1]])
#'
curve_ball<-function(g,nsim=1000,istrength=FALSE){
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

  if(istrength) {
    m <- get.adjacency(g,sparse=FALSE,attr="weight")
    r <- dim(m)[1]
    c <- dim(m)[2]

    wnets<- lapply(nets, function (x) {
      m1 <- get.adjacency(x,sparse=FALSE)
      for(i in seq_len(c) ){
        ss <- sample(m[,i])
        ss <- ss[ss>0]
        k <- 1
        for( j in seq_len(r)){
          if(m1[j,i]>0 ) {
            m1[j,i] <- ss[k]
            k <- k+1
          }
        }
      }
      graph_from_adjacency_matrix(m1,mode="directed",weighted = "weight")
    })
  } else {
    return(nets)
  }
}

#' @export
curveBall<-function(g,nsim=1000,istrength=FALSE){curve_ball(g,nsim,istrength)}



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
#' @aliases generateERbasal
#'
#' @return a list with igraph objects
#' @export
#'
#' @examples
#'
#' generateERbasal(netData[[1]])
generate_er_basal <- function(ig,nsim=1000){
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

#' @export
generateERbasal <- function(ig,nsim=1000){generate_er_basal(ig,nsim)}


#' Generate Niche Model Food Web(s)
#'
#' This function generates one or multiple food webs using the **Niche Model** proposed by Williams & Martinez (2000).
#' The model assumes that each species has a niche value and consumes resources within a defined range.
#'
#' @param S Integer. The number of species in the community. Must be `S > 1`.
#' @param C Numeric. The connectance (fraction of realized links). Must be `0 < C ≤ 1`.
#' @param nsim Integer. The number of networks to generate (`nsim ≥ 1`).
#'
#' @return If `nsim = 1`, returns a **binary adjacency matrix (`S × S`)**.
#' If `nsim > 1`, returns a **list of `igraph` objects**, each representing a food web.
#' @export
#'
#' @references
#' Williams, R. J., and N. D. Martinez. 2000. Simple rules yield complex food webs. *Nature* 404:180–183.
#'
#' @examples
#' generate_niche(20, 0.1)        # Single adjacency matrix
#' generate_niche(20, 0.1, nsim=5) # List of 5 food webs as igraph objects
generate_niche <- function(S, C, nsim = 1) {
  # Validate inputs
  if (!is.numeric(S) || S <= 1 || S %% 1 != 0) stop("S must be an integer greater than 1.")
  if (!is.numeric(C) || C <= 0 || C > 1) stop("C must be a number between 0 and 1.")
  if (!is.numeric(nsim) || nsim < 1 || nsim %% 1 != 0) stop("nsim must be a positive integer.")

  generate_single_network <- function() {
    connected <- FALSE
    while (!connected) {
      # Assign niche values
      n.i <- sort(runif(S))  # Sorted niche values

      # Determine feeding range
      r.i <- rbeta(S, 1, ((1 / (2 * C)) - 1)) * n.i  # Feeding range
      c.i <- runif(S, r.i / 2, n.i)  # Feeding center

      # Initialize adjacency matrix
      a <- matrix(0, nrow = S, ncol = S)

      # Assign feeding relationships
      for (i in 2:S) {  # Skip basal species
        for (j in 1:S) {
          if (n.i[j] > (c.i[i] - 0.5 * r.i[i]) && n.i[j] < (c.i[i] + 0.5 * r.i[i])) {
            a[j, i] <- 1
          }
        }
      }

      # Check connectivity
      g <- igraph::graph_from_adjacency_matrix(a, mode = "directed")
      connected <- igraph::is_connected(g)
    }
    return(g)  # Return igraph
  }

  if (nsim == 1) {
    return(generate_single_network())  # Return a single adjacency matrix
  } else {
    return(lapply(seq_len(nsim), function(x) generate_single_network()))  # Return list of igraph objects
  }
}
