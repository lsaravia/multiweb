
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


#' Potter Cove Body Mass and Biomass Dataset
#'
#' This dataset contains curated body mass, biomass, and resource density information
#' for consumers and resources in the Potter Cove food web, from Rodriguez & Saravia (2024).
#' The dataset was processed from raw measurements, with units converted to kilograms (kg),
#' and some special handling for detritus and sponge taxa where direct body mass or density values
#' were not available or meaningful.
#'
#' When resource density is unknown, the value `-999` is used as a flag so that
#' downstream functions can infer density or body mass using allometric relationships
#' following Pawar et al. (2012).
#' Body mass and biomass were originally in grams and converted to kilograms.
#' Sponge taxa (e.g., *Chalinidae*, *Dendrilla antarctica*) resource body mass and density
#' were assigned to 1.
#'
#' @format A data frame with the following columns:
#' \describe{
#'   \item{consumer}{Name of the consumer species.}
#'   \item{resource}{Name of the resource species.}
#'   \item{c_bodymass}{Body mass of the consumer (in kg).}
#'   \item{r_bodymass}{Body mass of the resource (in kg). May be set to `1` for detritus or sponges.}
#'   \item{r_biomass}{Biomass of the resource (in kg).}
#'   \item{r_density}{Resource density (units depend on context). May be `1` for some resources or `-999` to flag missing values.}
#' }
#'
#' @source <https://doi.org/10.5281/zenodo.10790590>
#'
#' @references
#' Rodriguez, I. D., & Saravia, L. A. (2024). Potter Cove’s Heavyweights: Estimation of Species’ Interaction Strength of an Antarctic Food Web. Ecology and Evolution, 14(11), e70389. https://doi.org/10.1002/ece3.70389
#'
#' Pawar, S., Dell, A.I., & Savage, V.M. (2012). Dimensionality of consumer search space drives trophic interaction strengths. \emph{Nature}, 486(7404), 485–489.
#'
#' @usage
#' PotterCove_bm
"PotterCove_bm"

