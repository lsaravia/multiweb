
#' Calculates the QSS difference between the full network and the network minus
#' one species
#'
#' The QSS determines the maximum eingenvalue of the community matrix (Jacobian) and characterizes
#' the linear stability of the network. This uses the function [multiweb::calc_QSS()] so it can take
#' into account the interaction strength if weigths are present. The comparison is made using the Anderson-Darling test with
#' the function [kSamples::ad.test()] and the Kolmogorov-Smirnov test [stats::ks.test()], both the p-values are reported as a
#' measure of strength of the diference.
#'
#' @param g        igraph network
#' @param sp_list  list with the species/nodes we will delete for the comparison
#' @param nsim     number of simulations to calculate QSS
#' @param ncores   number of cores used to perform the operation
#' @param istrength if TRUE takes the weigth attribute of the network as interaction strength to
#'                  calculate QSS.
#'
#' @return a data.frame with:
#'   * the node deleted
#'   * the Anderson-Darling statistic
#'   * the Anderson-Darling p-value
#'   * Kolmogorov-Smirnov statistic
#'   * Kolmogorov-Smirnov p-value
#'   * median of QSS of the complete network
#'   * median of QSS of the network with the deleted node
#'   * difference between the two QSS
#'
#' @import igraph
#' @importFrom dplyr bind_rows
#' @export
#'
#' @examples
#' \dontrun{
#' g <- netData[[1]]
#'
#' # Generate random weights
#' #
#' V(g)$weight <-  runif(vcount(g))
#'
#' # Without interaction strength
#' #
#' calc_QSS_extinction_dif(g,V(g)$name[1:3],nsim=10,istrength = FALSE)
#'
#' # With interaction strength
#' #
#' calc_QSS_extinction_dif(g,V(g)$name[1:3],nsim=10,istrength = TRUE)
#'}

calc_QSS_extinction_dif <- function(g, sp_list,nsim=1000, ncores=4, istrength = FALSE){

  # QSS for complete and deleted network
  QSS_all <- multiweb::calc_QSS(g, nsim = nsim, ncores = ncores, istrength = istrength, returnRaw = TRUE) %>%
    mutate(Network = "All spp")

  comp_webs <- lapply(sp_list, function(i){
    # delete one sp and create igraph object
    g_ext <- delete_vertices(g, i)
    size <- vcount(g_ext)
    # subset mean strength for deleted sp
    QSS_ext <- multiweb::calc_QSS(g_ext, nsim = nsim, ncores = ncores, istrength = istrength, returnRaw = TRUE) %>%
      mutate(Network = "One sp less")
    QSS <- bind_rows(QSS_all, QSS_ext)
    # extract p-value for Anderson-Darling test comparing complete and deleted network
    ad_test <- kSamples::ad.test(maxre ~ Network, data = QSS)
    ks_test <- ks.test(QSS_all$maxre,QSS_ext$maxre)
    # data frame
    data.frame(Deleted = i, Ad_stat = ad_test$ad[1,2], Ad_pvalue=ad_test$ad[1,3],
               KS_stat=ks_test$statistic,
               KS_pvalue=ks_test$p.value, QSS_all=median(QSS_all$maxre), QSS_ext=median(QSS_ext$maxre), difQSS = median(QSS_all$maxre)-median(QSS_ext$maxre))

  })
  bind_rows(comp_webs)
}




#' Calculates the QSS for an extinction sequence
#'
#' This functions calculates the QSS [multiweb::calc_QSS()] for a sequence of incremental deletions of species (nodes)
#' given by the vector `seq`, it does not produce secundary extinctions
#'
#' @param g_del igraph object for the deletion sequence
#' @param seq   vector or list with the nodes for the extinction sequence, the
#'              extinctions are incremental.
#' @param nsim        number of simulations used in QSS calculation
#' @param ncores      number of cores used to calculate QSS
#' @param istrength   parameter istrength for QSS function
#'
#' @return data frame with the number of remaining nodes, the connectance, the number unconnected of components, the median QSS, and the name
#'         of the last deleted node. The first row are the values for the complete network.
#' @export
#' @import igraph
#' @importFrom dplyr %>% summarize mutate
#'
#' @examples
#' \dontrun{
#' g <- netData[[1]]
#'
#' # Generate random weights
#' #
#' V(g)$weight <-  runif(vcount(g))
#'
#' # Without interaction strength
#' #
#' calc_QSS_extinctions_seq(g,V(g)$name[1:3],nsim=10,istrength = FALSE)
#'
#' # With interaction strength
#' #
#' calc_QSS_extinctions_seq(g,V(g)$name[1:3],nsim=10,istrength = TRUE)
#'}

calc_QSS_extinctions_seq <- function(g_del, seq, nsim=100, ncores=0, istrength=FALSE){
  lseq <- length(seq)
  if(lseq ==vcount(g_del)){                   # if the secuence is equal to the number of nodes delete the last component
    seq <- seq[1:(lseq-1)]
  }
  con <- multiweb::calc_topological_indices(g_del)$Connectance
  comp <- multiweb::calc_topological_indices(g_del)$Components
  totalqss <- multiweb::calc_QSS(g_del, nsim = nsim, ncores = ncores, istrength = istrength, returnRaw = TRUE) %>% summarize(QSS_median = median(maxre), QSS_prop = sum(maxre <0)/nsim) %>%
    mutate(Size = vcount(g_del), Connectance = con, Components = comp, Last_deleted = "Total")

  qdel <- lapply(seq, function(i){
    g_del <<- delete_vertices(g_del, i)
    size <- vcount(g_del)
    con <- multiweb::calc_topological_indices(g_del)$Connectance
    comp <- multiweb::calc_topological_indices(g_del)$Components
    QSS <- multiweb::calc_QSS(g_del, nsim = nsim, ncores = ncores, istrength = istrength, returnRaw = TRUE) %>% summarize(QSS_median = median(maxre), QSS_prop = sum(maxre <0)/nsim) %>%
      mutate(Size = size, Connectance = con, Components = comp, Last_deleted = i)
  })
  qdel <- bind_rows(totalqss,qdel)
}

