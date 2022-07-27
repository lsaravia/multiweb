#' Calcultates the interaction intensity of a food web using the metabolic theory and
#' the interaction dimensionality
#'
#' It uses all the coeficients from [1], taking into account the interaction dimensionality and body mass units in Kg.
#' For detritus or sediment the resource mass mean `res_mm` is generally < 0 thus resource body mass is calculated using
#' the equation S9 and supplementary figures 2c & d (kg) of the paper.
#' The values for the resource density $Xr$ where estimated according the equation S18 and supplementary figures 2i & j (individuals/m2 - m3)
#'
#' @references
#' 1. Pawar, S., Dell, A. I., & Van M. Savage. (2012). Dimensionality of consumer search space drives trophic interaction strengths. Nature, 486, 485. https://doi.org/10.1038/nature11131
#'
#' @param da data.frame with the interactions body mass and type of interaction dimensionality
#' @param res_mm name of the column with the resource mass mean
#' @param con_mm name of the column with the consumer mass mean
#' @param int_dim name of the column with the interaction dimensionality
#'
#' @return a data.frame
#'
#' @import dplyr
#'
#' @export
#'
#' @examples
calc_interaction_intensity <- function(da,res_mm,con_mm,int_dim) # alfa0 = alfa2D/3D Pa = exponente of alfa 2D/3D
{
  #  con_taxonomy            res_taxonomy              con_mass_mean res_mass_mean interaction_dim
  #
  #
  # if( typeof(wedd_df$interaction_dim) != "character")
  #   stop("int_dim must be character type")
  #
  det <- da %>% filter(is.na({{int_dim}}))
  if( nrow(det) > 0 ) warning(paste("Interaction dimensionality is not defined for", nrow(det), "rows"))
  d2D <- filter(da, {{int_dim}}=="2D") %>%
    mutate(mR=if_else({{res_mm}} < 0, 10^-2.6*{{con_mm}}^0.67 , {{res_mm}} ),
           xR= 10^-2.67 * mR^-0.79,
           alfa=10^-3.08 *{{con_mm}}^0.68,
           qRC = alfa*xR*mR/{{con_mm}})

  d3D <- filter(da, {{int_dim}}=="3D") %>%
    mutate(mR=if_else({{res_mm}} < 0, 10^-2.96*{{con_mm}}^1.46 , {{res_mm}} ),
           xR= 10^-2.48 * mR^-0.86,
           alfa = 10^-1.77 * {{con_mm}}^1.05,
           qRC = alfa*xR*mR/{{con_mm}}
    )
  bind_rows(d2D,d3D)

}
