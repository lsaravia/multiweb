
#' Metadata for `netData`
#'
#' A data.frame with extended information about the 29 networks from (1), included in `netData`
#'
#' @format A data frame
#' \describe{
#'   \item{Network}{Name of the network}
#'   \item{Longitude}{Longitude of the geographical coordinates of the network}
#'   \item{Latitude}{Latitude of the geographical coordinates of the network}
#'   \item{Date}{Date range or approximate date when the data was collected}
#'   \item{InteractionsSource}{Method used to determine the interactions}
#'   \item{Reference}{Reference of the study from where data was obtained}
#'   \item{ReferenceLink}{URL of the reference}
#' }
#'
#' @source \url{https://doi.org/10.1371/journal.pone.0198217}
#'
#' 1. Marina, T. I., Saravia, L. A., Cordone, G., Salinas, V., Doyle, S. R., & Momo, F. R. (2018). Architecture of marine food webs: To be or not be a ‘small-world.’ PLoS ONE, 13(5), 1–13. https://doi.org/10.1371/journal.pone.0198217
#'
"metadata"


#' Compilation of Ecological Networks (food webs) in igraph format.
#'
#' A dataset 29 curated and highly-resolved marine networks (food webs), updated from (1)
#'
#' @format A list with 29 igraph network objects
#'
#' @source \url{https://doi.org/10.1371/journal.pone.0198217}
#'
#' @references
#'
#' 1. Marina, T. I. et al. 2018. Architecture of marine food webs: To be or not be a ‘small-world.’ - PLoS ONE 13: 1–13.
"netData"
