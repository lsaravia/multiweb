#' Calcultates the interaction intensity of a food web using the metabolic theory and
#' the interaction dimensionality
#'
#' It uses the body mass in Kg of predator/consumer and prey/resources and the dimensionality of the interaction as source data,
#' then with all the coeficients from ref. 1.
#' For detritus or sediment the resource mass mean `res_mm` is generally < 0 thus resource body mass (kg) is calculated using
#' the equation S9 and supplementary figures 2c & d of the paper.
#' The values for the resource density *Xr* where estimated according the equation S18 and supplementary figures 2i & j (individuals/m2 - m3)
#'
#' @references
#' 1. Pawar, S., Dell, A. I., & Van M. Savage. (2012). Dimensionality of consumer search space drives trophic interaction strengths. Nature, 486, 485. https://doi.org/10.1038/nature11131
#'
#' @param da data.frame with the interactions body mass and type of interaction dimensionality
#' @param res_mm name of the column with the resource mass mean
#' @param con_mm name of the column with the consumer mass mean
#' @param int_dim name of the column with the interaction dimensionality
#'
#' @return A data.frame based on `da` with the following fields added
#'
#' * mR:  if `res_mm<0` is the resource mass calculated with the equations from ref 1,
#'        if `res_mm>0` duplicates the value of `res_mm`
#' * xR: calculated resource density.
#' * alfa: calculated search rate.
#' * qRC: calculated trophic interaction strength as `alfa*xR*mR/mC` where `mC` is the consumer body mass.
#'
#' @importFrom dplyr %>% filter mutate if_else bind_rows
#' @export
#'
#' @examples
#' \dontrun{
#' g <- netData[[1]]
#'
#' require(dplyr)
#'
#' # build the data.frame with random values
#'
#' da <- as_long_data_frame(g) %>% dplyr::select(from:to) %>% mutate(con_mm=rlnorm(n(),5,2),res_mm=con_mm - 30 ,int_dim=sample(c("2D","3D"),n(),replace=TRUE))
#'
#' calc_interaction_intensity(da,res_mm,con_mm,int_dim)
#' }

calc_interaction_intensity <- function(da,res_mm,con_mm,int_dim) # alfa0 = alfa2D/3D Pa = exponente of alfa 2D/3D
{
  #  con_taxonomy            res_taxonomy              con_mass_mean res_mass_mean interaction_dim
  #
  #
  # if( typeof(wedd_df$interaction_dim) != "character")
  #   stop("int_dim must be character type")
  #
  # resource Mass vs consumer Mass Exponents (Sup Fig 2 )
  pm2D = 0.67 # sd = 0.17
  pm3D = 1.46 # sd = 0.12
  int_pm2D = -2.6   # intercepts
  int_pm3D = -2.96
  #
  det <- da %>% filter(is.na({{int_dim}}))
  if( nrow(det) > 0 ) warning(paste("Interaction dimensionality is not defined for", nrow(det), "rows"))
  d2D <- filter(da, {{int_dim}}=="2D") %>%
    mutate(mR=if_else({{res_mm}} < 0, 10^int_pm2D*{{con_mm}}^pm2D , {{res_mm}} ),
           xR= 10^-2.67 * mR^-0.79,
           alfa=10^-3.08 *{{con_mm}}^0.68,
           qRC = alfa*xR*mR/{{con_mm}})

  d3D <- filter(da, {{int_dim}}=="3D") %>%
    mutate(mR=if_else({{res_mm}} < 0, 10^int_pm3D*{{con_mm}}^pm3D , {{res_mm}} ),
           xR= 10^-2.48 * mR^-0.86,
           alfa = 10^-1.77 * {{con_mm}}^1.05,
           qRC = alfa*xR*mR/{{con_mm}}
    )
  bind_rows(d2D,d3D)
#
}
