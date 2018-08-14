
#' Curbe ball algorithm to generate random networks from a given network with the same colums and row totals
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

  nets<- lapply(1:nsim, function (x) {

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

    graph_from_adjacency_matrix(rm,mode="directed")

  })
}
