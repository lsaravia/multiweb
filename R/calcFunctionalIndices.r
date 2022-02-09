#' Calculates the quantitative connectance, and the effective number of flows using the interaction matrix and a vector of species abundances
#'
#' The Quantitative connectance (Cq) takes into account both the distribution of per capita interaction strengths
#' among species in the web and the distribution of species’ abundances and quantifies the diversity of network fluxes,
#' if all the species have the same flux is equal to the directed connectance.
#' The mean or effective number of flows impinging upon or emanating from a tipical node (LDq) is based the average flow diversity, and when all flows are equal is similar to linkage density. Both measures are based in Shannon information theory.
#' The total interaction flux is measured
#' as ```T[i,j] <- d[i] * d[j] * interM[i,j]```. The effective Cq is calculated following the formulas in appendix 2 of [1], LDq follows [2]
#'
#' @references
#' 1. Fahimipour, A.K. & Hein, A.M. (2014). The dynamics of assembling food webs. Ecol. Lett., 17, 606–613
#'
#' 1. Ulanowicz, R.E. & Wolff, W.F. (1991). Ecosystem flow networks: Loaded dice? Math. Biosci., 103, 45–68
#'
#' @param interM per capita interaction strength matrix
#' @param d      species' abundances vector
#'
#' @return A list with Cq,the quantitative connectance index and LDq, =
#'
#' @aliases calcQuantitativeConnectance
#'
#' @export
#'
#' @examples
#'
#' # 3 predators 2 preys unequal fluxes
#' #
#' m <- matrix()
#' matrix(0,nrow=4,ncol=4)
#' m[1,2] <- m[1,3] <- m[3,4]<- .2
#' m[2,1] <- m[3,1] <- m[4,3] <- -2
#'
#' calc_quantitative_connectance(m, c(1,1,1,1))
#'
#'# Equal input and output fluxes
#'
#'m <- matrix(0,nrow=4,ncol=4)
#'m[1,2] <- m[1,3] <- m[3,4]<- 2
#'m[2,1] <- m[3,1] <- m[4,3] <- -2
#'calc_quantitative_connectance(m, c(1,1,1,1))


calc_quantitative_connectance <- function(interM,d){

  nsp <- length(d) # number of species
  if(ncol(interM)!=nrow(interM)) stop("interM has to be square, number of rows has to be equal to number of columns")
  if(nsp!=nrow(interM)) stop("length of species vector has to be equal to interM dimensions")

  TT <- matrix(0,nsp,nsp)
  SS <- matrix(0,nsp,nsp)

  Hout <- numeric(nsp)
  Hin  <- numeric(nsp)

  for(i in seq_along(d)) {
    for(j in seq_along(d)) {
      TT[i,j] <- abs(interM[i,j])*d[i]*d[j]
    }
  }
  TTcol <- colSums(TT)    # Input
  TTrow <- rowSums(TT)    # Output
  TTT <- sum(TTrow)

  for(k in seq_along(d)) {
    for(i in seq_along(d)) {
      SS[i,k] <- TT[i,k]/TTT * log2(TT[i,k]^2/(TTrow[i]*TTcol[k]))
      if(TT[i,k]!=0){
        Hin[k] <- Hin[k] - TT[i,k]/TTcol[k]*log(TT[i,k]/TTcol[k])
      }
      if(TT[k,i]!=0){
        Hout[k] <- Hout[k] - TT[k,i]/TTrow[k]*log(TT[k,i]/TTrow[k])
      }
    }
  }
  # SS1 is the \Sigma separated in input and output components see Ulanowicz eq 3
  # (SS1 <- (sum(TTcol*Hin)+sum(TTrow*Hout))/TTT) == sum(SS)
  #
  nout <- exp(Hout)
  nout[TTrow==0] <- 0
  nin <- exp(Hout)
  nin[TTcol==0] <- 0
  LDq <- 2^(-sum(SS,na.rm = T)/2)   # C= linkage Density / S
  Cq <- 1/(2*nsp^2)*(sum(nout+nin))  #
  return(list(Cq=Cq,LDq=LDq))
}

calcQuantitativeConnectance <- function(interM,d){
  calc_quantitative_connectance(interM,d)}


#' Calc the Quasi Sign Stability measure for antagonistic (predator-prey) or mgraph networks with multiple interactions.
#'
#' The proportion of matrices that are locally stable, these matrices are created by sampling the values of the community matrix
#' (the Jacobian) from a uniform distribution, preserving the sign structure [1]. If the 'ig' parameter is
#' an `mgraph` network it needs to have been built with the order `c("Competitive", "Mutualistic", "Trophic")`
#' It also calculates the mean of the real part of the maximum eingenvalue, which is also a measure of stability [2].
#' It uses a uniform distribution between 0 and maximum values given by the parameters `negative`, `positive` and `selfDamping`,
#' corresponding to the sign of interactions and self-limitation effect [3,4].
#' If the edges of the networks have a weigth attribute and `istrength` parameter is true, weigth will be used as interaction strength, then the limits of the uniform distribution
#' will be `negative*-x`, `positive*x` where x is the value of the weigth for the edge.
#' If the values of these parameters are 0 then there is no interaction of that kind.
#'
#' @references
#'
#' 1. Allesina, S. & Pascual, M. (2008). Network structure, predator - Prey modules, and stability in large food webs.
#' Theor. Ecol., 1, 55–64.
#' 2. Grilli, J., Rogers, T. & Allesina, S. (2016). Modularity and stability in ecological communities. Nat. Commun., 7, 12031
#' 3. Monteiro, A.B. & Del Bianco Faria, L. (2017). Causal relationships between population stability and food-web topology. Functional Ecology, 31, 1294–1300.
#' 4. Borrelli, J. J. 2015. Selection against instability: stable subgraphs are most frequent in empirical food webs. - Oikos 124: 1583–1588.


#'
#' @param ig  igraph or a list of igraph networks or mgraph network
#' @param nsim number of simulations to calculate QSS
#' @param ncores number of cores to use in parallel comutation if 0 it uses sequential processing
#' @param negative the maximum magnitude of the negative interaction (the effect of the predator on the prey) must be <= 0
#' @param positive the maximum magnitude of the positive interaction (the effect of the prey on the predator) must be >= 0
#' @param selfDamping the maximum magnitude of the self-limitation (the effect of the species on itself) must be <= 0
#' @param istrength If TRUE takes the weigth attribute of the network as interaction strength.
#'
#' @return a data.frame with the QSS, and MEing, the mean of the real part of the maximum eingenvalue
#'
#' @export
#'
#' @import igraph
#' @importFrom future sequential multiprocess
#' @importFrom future.apply future_lapply
#' @importFrom dplyr bind_rows
#'
#' @examples
#' \dontrun{
#'
#' g <- netData[[2]]
#'
#' tp <- calc_QSS(g)
#'
#' # Read Multiplex network and calculate QSS
#'
#' fileName <- c(system.file("extdata",  package = "multiweb"))
#' dn <- list.files(fileName,pattern = "^Kefi2015.*\\.txt$")
#' gt <- readMultiplex(dn,types,"inst/extdata", skipColumn = 2)
#' calc_QSS(gt)
#'
#' }
calc_QSS <- function(ig,nsim=1000,ncores=0,negative=-10, positive=0.1, selfDamping=-1,istrength=FALSE) {
  if(inherits(ig,"igraph")) {
    ig <- list(ig)
  } else if(class(ig[[1]])!="igraph") {
    stop("parameter ig must be an igraph object")
  }
  if(negative>0 || positive<0 || selfDamping>0)
    stop("Parameters should be negative<=0, positive>=0, selfDamping<=0")

  if(ncores) {
    cn <- future::availableCores()
    if(ncores>cn)
      ncores <- cn
    future::plan(multiprocess, workers=ncores)
    on.exit(future::plan(sequential))
  } else {
    future::plan(sequential)
  }

  if(class(ig)!='mgraph') {
    df <-  lapply(ig, function(red)
    {
      lred <- fromIgraphToMgraph(list(make_empty_graph(n=vcount(red)),make_empty_graph(n=vcount(red)),red),
                                   c("empty","empty","Antagonistic"))
      mat <- toGLVadjMat(lred,c("empty","empty","Antagonistic"),istrength = istrength)   #

      df <- future_lapply(seq_len(nsim), function(i){
        ranmat <- ranUnif(mat,negative,positive,selfDamping)
        eigs <- maxRE(ranmat)
      },future.seed = TRUE)
      df <-   do.call(rbind,df)
      data.frame(QSS=sum(df<0)/nsim,MEing=mean(df))
    })
  } else {
    mat <- toGLVadjMat(ig,c("Competitive", "Mutualistic", "Trophic"),istrength = istrength)   #
    df <- future_lapply(seq_len(nsim), function(i){
      ranmat <- ranUnif(mat,negative,positive,selfDamping)
      eigs <- maxRE(ranmat)
    },future.seed = TRUE)
    df <-   do.call(rbind,df)
    df <- data.frame(QSS=sum(df<0)/nsim,MEing=mean(df))
  }

  bind_rows(df)
}


# Auxiliar function of calc_QSS, Calculates a random community mattrix with a fixed signed structure
#
ranUnif <- function(motmat, negative=-10,positive=0.1,selfDamping=-1){
  newmat <- apply(motmat, c(1,2), function(x){
    #if(x==1){runif(1, 0, positive)}else if(x==-1){runif(1, negative, 0)} else{0}
    if(x>0){runif(1, 0, x*positive)}else if(x<0){runif(1,-x*negative, 0)} else{0}
  })
  diag(newmat) <- runif(nrow(motmat), selfDamping, 0)
  return(newmat)
}

# Auxiliar function of calc_QSS, Compute the eigenvalues and return the largest real part.
#
maxRE <- function(rmat){
  lam.max <- max(Re(eigen(rmat)$values))
  return(lam.max)
}


#' Function to calculate weighted network indices
#'
#' The function calculates: weighted linkage density, connectance, generality and vulnerability and the SDs of the last two based on [1]
#' and the level of omnivory, mean and maximum trophic level based on [2] using [NetIndices::TrophInd()] function.
#' The igraph networks must have the weight attribute.
#'
#' @param ig An igraph object with the weight attribute representing fluxes
#' @param ncores number of cores used to compute in parallel, if 0 sequential processing is used.
#'
#' @references
#'
#' 1. Bersier, LF. et al. (2002). Quantitative Descriptors Of Food-web Matrices. Ecology, 83(9), 2394–2407.
#' 2. Kones, JK. et al. (2009). Are network indices robust indicators of food web functioning? A Monte Carlo approach.
#'    Ecological Modelling, 220, 370–382.
#'
#'
#' @return a data.frame with the following fields:
#'
#'  \item{LD:}{ linkage density}
#'  \item{Connectance:}{ directed Connectance}
#'  \item{TLmean:}{ mean trophic level}
#'  \item{TLmax:}{ maximum trophic level}
#'  \item{LOmnivory:} { Level of omnivory, quantiﬁes mean of the variety in trophic levels of the preys of a consumer}
#'  \item{Vulnerability:}{ mean of number of consumers per prey}
#'  \item{VulSD:}{ the standard deviation of normalized Vulnerability}
#'  \item{Generality:}{ mean number of prey per consumer}
#'  \item{GenSD:}{ the standard deviation of normalized Generality}
#'
#'
#' @export
#' @importFrom NetIndices TrophInd
#' @import igraph
#' @importFrom foreach foreach %dopar%
#' @importFrom doFuture registerDoFuture
#' @importFrom future sequential multiprocess
#'
#' @examples
#'
#' # Generate a test network
#'
#' g <-   graph_from_literal( 1 -+ 4 -+ 7,2 -+ 5 -+7, 3-+6-+7, 7-+7, 4+-3, simplify = FALSE)
#'
#' # Add weight
#'
#' E(g)$weight <- sample(c(.1,.2,.8,.9),gsize(g),replace=TRUE)
#'
#' calc_weighted_topological_indices(g)

calc_weighted_topological_indices<- function(ig,ncores=0){

  if(inherits(ig,"igraph")) {
    ig <- list(ig)
  } else if(class(ig[[1]])!="igraph") {
    stop("parameter ig must be an igraph object")
  }

  registerDoFuture()
  if(ncores) {
    cn <- future::availableCores()
    if(ncores>cn)
      ncores <- cn
    future::plan(multiprocess, workers=ncores)
    on.exit(future::plan(sequential))
  } else {
    future::plan(sequential)
  }


  df <-  foreach(g=ig,.combine='rbind',.inorder=FALSE,.packages=c('igraph','NetIndices')) %dopar% {

    W.net <- as_adjacency_matrix(g,sparse=FALSE,attr="weight")
    res<-c()
    # The flux matrix

    ### Taxon-specific Shannon indices of inflows
    # sum of k species inflows --> colsums
    sum.in<-apply(W.net, 2, sum)

    # Diversity of k species inflows
    H.in.mat<- t(t(W.net)/sum.in)*t(log(t(W.net)/sum.in)) #columns divided by the total col sum
    H.in.mat[!is.finite(H.in.mat)] <- 0 #converts NaN to 0's
    H.in<- apply(H.in.mat, 2, sum)*-1

    ### Taxon-specific Shannon indices of outflows
    # sum of k speies outflows --> rowsums
    sum.out<-apply(W.net, 1, sum)

    # Diversity of k species outflows
    H.out.mat<- (W.net/sum.out)*log(W.net/sum.out) #rows divided by the total row sum
    H.out.mat[!is.finite(H.out.mat)] <- 0 #converts NaN to 0's
    H.out<- apply(H.out.mat, 1, sum)*-1

    # Effective number of prey or resources = N(R,k)
    # The reciprocal of H(R,k) --> N (R,k) is the equivalent number of prey for species k
    # N.res<-exp(H.in)
    N.res<-ifelse(sum.in==0, H.in, exp(H.in))

    # Effective number of predators or consumers = N(C,k)
    # The reciprocal of H(C,k) --> N (C,k) is the equivalent number of predators for species k
    # N.con<-exp(H.out)
    N.con<-ifelse(sum.out==0, H.out, exp(H.out))

    ### Quantitative Weighted and unweighted link density and weighted connectance
    no.species<-ncol(W.net)

    # The weighted link density (LDw) is:
    # In the weighted version the effective number of predators for species i is weighted by i's
    # contribution to the total outflow
    # the same is the case for the inflows

    tot.mat<- sum(W.net)
    # LD.w <- (sum((sum.in/tot.mat)*N.res) + sum((sum.out/tot.mat)*N.con))/2
    # equivalent to next formula, but next one is closer to manuscript
    res$qLD.w <- 1/(2*tot.mat)*(sum(sum.in*N.res) + sum(sum.out*N.con))

    # Weighted connectance
    res$qC.w<- res$qLD.w/no.species

    # defintion according to Bersier et al. 2002 top species = [0.99, 1]
    # int.sp<- pos.ind[pos.ind>0&pos.ind<0.99]# intermediate are all species that are not basal nor top

    # weighted quantitative Generality
    res$qG.w<-sum(sum.in*N.res/sum(W.net))

    # weighted quantitative Vulnerability
    res$qV.w<-sum(sum.out*N.con/sum(W.net))

    # Standard deviation of stanardized unweighted quantitative Generality
    s<-no.species # length(pos.ind)

    # Standard deviation of stanardized weighted quantitative Generality
    res$qGsd.w<-sd(s*sum.in*N.res/sum(N.res*sum.in))

    # Standard deviation of stanardized weighted quantitative Vulnerability
    res$qVsd.w<-sd(s*sum.out*N.con/sum(N.con*sum.out))

    # Weigthed trophic levels
    #
    wtl <- TrophInd(W.net)

    wOmn <- mean(wtl$OI)        # Omnivory Level

    data.frame( wLD=res$qLD.w, wConnectance=res$qC.w, wTLmean=mean(wtl$TL),wTLmax=max(wtl$TL),wLOmnivory=wOmn,
                wVulnerability=res$qV.w,wVulSD=res$qVsd.w,wGenerality=res$qG.w,wGenSD=res$qGsd.w)
  }
  return(df)

}
