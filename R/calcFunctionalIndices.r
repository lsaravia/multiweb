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
