
#' Calculates the QSS difference between the full network and the network minus
#' one species
#'
#' The Quasi-sign stability is estimated with of the community matrix (Jacobian) and characterizes
#' the linear stability of the network. This uses the function [multiweb::calc_QSS()] so it can take
#' into account the interaction strength if weights are present. The comparison is made using the Anderson-Darling test with
#' the function [kSamples::ad.test()] and the Kolmogorov-Smirnov test [stats::ks.test()], both the p-values are reported as a
#' measure of strength of the difference. If istrength is TRUE it makes a comparison to a null model with the same species and links than the
#' reduced network but with all interaction strengths equal to the mean interaction strength.
#'
#' @param g        igraph network
#' @param sp_list  list with the species/nodes we will delete for the comparison
#' @param nsim     number of simulations to calculate QSS
#' @param ncores   number of cores used to perform the operation
#' @param istrength if TRUE takes the weight attribute of the network as interaction strength to
#'                  calculate QSS.
#'
#' @return a data.frame with:
#'   * the node deleted
#'   * the Anderson-Darling p-value of the comparison with the complete network
#'   * Kolmogorov-Smirnov p-value of the comparison with the complete network
#'   * If `istrength == TRUE` the Anderson-Darling p-value of the comparison with the null network
#'   * If `istrength == TRUE` Kolmogorov-Smirnov p-value of the comparison with the null network
#'   * median of QSS of the complete network
#'   * median of QSS of the network with the deleted node
#'   * difference between the two previous median QSS
#'   * If `istrength == TRUE` median of QSS of the null network
#'   * If `istrength == TRUE` difference between the median QSS of the network with the deleted node
#'     and the null network.
#'
#'
#' @import igraph
#' @importFrom dplyr bind_rows
#' @importFrom kSamples ad.test
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

  if( istrength ) {
    if( length(edge_attr_names(g))==0 ) {
      stop("Network must have weight attribute")
    } else {
      if( !any(edge_attr_names(g) %in% "weight") )
        stop("Network must have weight attribute")
    }
  }

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

    if( istrength ) {
      # Build Null model with equal mean interaction strength
      #
      mean_w <- mean(E(g_ext)$weight)
      g1 <- g_ext
      E(g1)$weight <- mean_w
      QSS_null <- multiweb::calc_QSS(g1, nsim = nsim, ncores = ncores, istrength = istrength, returnRaw = TRUE) %>%
        mutate(Network = "One sp less null")
      QSS <- bind_rows(QSS_ext, QSS_null)

      ad_testn <- kSamples::ad.test(maxre ~ Network, data = QSS)
      ks_testn <- ks.test(QSS_ext$maxre,QSS_null$maxre)

      data.frame(Deleted = i,
                 Ad_pvalue=ad_test$ad[1,3],
                 KS_pvalue=ks_test$p.value,
                 Ad_pvalue_null=ad_testn$ad[1,3],
                 KS_pvalue_null=ks_testn$p.value,
                 QSS_all=median(QSS_all$maxre), QSS_ext=median(QSS_ext$maxre), difQSS = median(QSS_all$maxre)-median(QSS_ext$maxre),
                 QSS_null=median(QSS_null$maxre) ,difQSS_null = median(QSS_ext$maxre)-median(QSS_null$maxre)
                 )

    } else {
    # data frame

      data.frame(Deleted = i,
               Ad_pvalue=ad_test$ad[1,3],
               KS_pvalue=ks_test$p.value,
               QSS_all=median(QSS_all$maxre), QSS_ext=median(QSS_ext$maxre), difQSS = median(QSS_all$maxre)-median(QSS_ext$maxre)
        )
    }
  })
  bind_rows(comp_webs)
}


#' Calculates the QSS difference between the full network and the network minus
#' a group of species
#'
#' The Quasi-sign stability is estimated with the maximum eigenvalue of the community matrix (Jacobian) and characterizes
#' the linear stability of the network. This uses the function [multiweb::calc_QSS()] so it can take
#' into account the interaction strength if weights are present. The comparison is made using the Anderson-Darling test with
#' the function [kSamples::ad.test()] and the Kolmogorov-Smirnov test [stats::ks.test()], both the p-values are reported as a
#' measure of strength of the difference.
#'
#' @param g        igraph network
#' @param sp_list  list with the group pf species/nodes we will delete for the comparison
#' @param nsim     number of simulations to calculate QSS
#' @param ncores   number of cores used to perform the operation
#' @param istrength if TRUE takes the weight attribute of the network as interaction strength to
#'                  calculate QSS.
#'
#' @return a data.frame with:
#'   * the Anderson-Darling p-value of the comparison with the complete network
#'   * Kolmogorov-Smirnov p-value of the comparison with the complete network
#'   * If `istrength == TRUE` the Anderson-Darling p-value of the comparison with the null network
#'   * If `istrength == TRUE` Kolmogorov-Smirnov p-value of the comparison with the null network
#'   * median of QSS of the complete network
#'   * median of QSS of the network with the deleted nodes
#'   * difference between the two previous median QSS
#'
#' @import igraph
#' @importFrom dplyr bind_rows
#' @importFrom kSamples ad.test
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
#' calc_QSS_extinction_dif_grp(g,V(g)$name[1:3],nsim=10,istrength = FALSE)
#'
#' # With interaction strength
#' #
#' calc_QSS_extinction_dif_grp(g,V(g)$name[1:3],nsim=10,istrength = TRUE)
#'}

calc_QSS_extinction_dif_grp <- function(g, sp_list,nsim=1000, ncores=4, istrength = FALSE){

  if( istrength ) {
    if( length(edge_attr_names(g))==0 ) {
      stop("Network must have weight attribute")
    } else {
      if( !any(edge_attr_names(g) %in% "weight") )
        stop("Network must have weight attribute")
    }
  }

  # QSS for complete and deleted network
  QSS_all <- multiweb::calc_QSS(g, nsim = nsim, ncores = ncores, istrength = istrength, returnRaw = TRUE) %>%
    mutate(Network = "All spp")

  # delete one sp and create igraph object
  g_ext <- delete_vertices(g, sp_list)
  con <- multiweb::calc_topological_indices(g_ext)$Connectance
  comp <- multiweb::calc_topological_indices(g_ext)$Components
  size <- vcount(g_ext)
  # Calculates the QSS
  QSS_ext <- multiweb::calc_QSS(g_ext, nsim = nsim, ncores = ncores, istrength = istrength, returnRaw = TRUE) %>%
    mutate(Network = "Minus Group")
  QSS <- bind_rows(QSS_all, QSS_ext)
  # extract p-value for Anderson-Darling test comparing complete and deleted network
  ad_test <- kSamples::ad.test(maxre ~ Network, data = QSS)
  ks_test <- ks.test(QSS_all$maxre,QSS_ext$maxre)

  # result data frame

  data.frame(Size = size, Connectance = con, Components = comp,
               Ad_pvalue=ad_test$ad[1,3],
               KS_pvalue=ks_test$p.value,
               QSS_all=median(QSS_all$maxre), QSS_ext=median(QSS_ext$maxre), difQSS = median(QSS_all$maxre)-median(QSS_ext$maxre)
    )

}

#' Calculates the QSS for an extinction sequence
#'
#' This functions calculates the QSS [multiweb::calc_QSS()] for a sequence of incremental deletions of species (nodes)
#' given by the vector `seq`, it does not produce secondary extinctions
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

