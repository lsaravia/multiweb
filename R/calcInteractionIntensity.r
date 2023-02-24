#' Calculates the interaction intensity of a food web using the metabolic theory and
#' the interaction dimensionality
#'
#' The function uses the body mass in Kg of predator/consumer and prey/resources and the dimensionality of the interaction as source data,
#' then the interacion intensity is estimated with all the coeficients from [1] as `alfa*xR*mR/mC`, where `alpha` is the search rate `xR`
#' the resource density, `mR` the resource body mass and `mC` the consumer body mass. This value of the interaction strength quantifies
#' the effect of the predator on the prey. Assuming a Lotka-Volterra model is equivalent to the entry A(i,j) of the Jacobian, where i is the
#' prey and j the predator.
#'
#' If the resource density is not known (parameter `res_den`) you must set the column to a less than 0 value; and it
#' will be estimated according to the equation S18 and supplementary figures 2i & j (individuals/m2 - m3)
#'
#' For detritus or sediment the resource mass mean is not known (parameter `res_mm`) it must be set as negative;
#' and the  resource body mass (kg) will be calculated using the equation S9 and supplementary figures 2c & d of the paper.
#'
#' If the parameter 'nsims > 1 ' the function will estimate the variability on each interaction strength. It takes random values from a normal distribution
#' with mean and standard deviation given by the Pawar's regressions for the slopes of allometric exponents.
#'
#' @references
#' 1. Pawar, S., Dell, A. I., & Van M. Savage. (2012). Dimensionality of consumer search space drives trophic interaction strengths. Nature, 486, 485. https://doi.org/10.1038/nature11131
#'
#' @param da data.frame with the interactions body mass and type of interaction dimensionality
#' @param res_mm name of the column with the resource body mass mean
#' @param res_den name of the column with the resource density in Individuals/m^2 in 2D or m^3 in 3D. If lower than 0 it uses the previously mentioned
#'                estimation.
#' @param con_mm name of the column with the consumer body mass mean
#' @param int_dim name of the column with the interaction dimensionality
#'
#' @return A data.frame based on `da` with the following fields added
#'
#' * mR:  if `res_mm<0` is the resource mass calculated with the equations from ref 1,
#'        if `res_mm>0` duplicates the value of `res_mm`
#' * xR: calculated resource density or the same value as in the input data.frame.
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
#' set.seed(7815)
#' da <- as_long_data_frame(g) %>% dplyr::select(from:to) %>% mutate(con_mm=rlnorm(n(),5,2),res_mm=con_mm - 30 ,int_dim=sample(c("2D","3D"),n(),replace=TRUE), res_den = -999)
#'
#' calc_interaction_intensity(da,res_mm,res_den,con_mm,int_dim, nsims=1)
#' }

calc_interaction_intensity <- function(da,res_mm,res_den,con_mm,int_dim, nsims=1)
{
  #
  # if( typeof(wedd_df$interaction_dim) != "character")
  #   stop("int_dim must be character type")
  #
  # resource Mass vs consumer Mass Exponents (Sup Fig 2 c&d)
  pm2D = 0.67 # 95% CI 0.17 sd = SE * sqrt(n) / 1.96, n = 39 (total n overestimate SD)
  pm3D = 1.46 # 95% CI 0.12
  int_pm2D = -2.6   # intercepts
  int_pm3D = -2.96
  # Minimun Resource Density
  #
  px2D <- -0.79 # 95% CI ± 0.09  Fig 2 i&j n = 123 (total n)
  px3D <- -0.86 # 95% CI ± 0.06            n = 131
  x02D <- -2.67
  x03D <- -2.48
  # Search Rate alfa scaling
  #
  al2D <- 0.68 # CI 0.12 n = 255 (we use the total n reported in the paper which overestimates the SD)
  al3D <- 1.05 # CI 0.08 n = 255
  int_al2D <- -3.08
  int_al3D <- -1.77

  det <- da %>% filter(is.na({{int_dim}}))
  if( nrow(det) > 0 ) warning(paste("Interaction dimensionality is not defined for", nrow(det), "rows"))

  if(nsims == 1 ){
  d2D <- filter(da, {{int_dim}}=="2D") %>%
    mutate(mR=if_else({{res_mm}} < 0, 10^int_pm2D*{{con_mm}}^pm2D , {{res_mm}} ),
           xR=if_else({{res_den}} < 0 , 10^x02D * mR^px2D,{{res_den}}),
           alfa=10^int_al2D *{{con_mm}}^al2D,
           qRC = alfa*xR*mR/{{con_mm}})

  d3D <- filter(da, {{int_dim}}=="3D") %>%
    mutate(mR=if_else({{res_mm}} < 0, 10^int_pm3D*{{con_mm}}^pm3D , {{res_mm}} ),
           xR=if_else({{res_den}} < 0 , 10^x03D * mR^px3D,{{res_den}}),
           alfa = 10^int_al3D * {{con_mm}}^al3D,
           qRC = alfa*xR*mR/{{con_mm}}
    )
  da <- bind_rows(d2D,d3D)
  } else {
    #
    #
    da <- lapply(seq_len(nsims), function(x){
      pm2D <- rnorm(1,0.67, 0.17*sqrt(19)/1.96)    # 95% CI 0.17 sd = SE * sqrt(n)/1.96, n = 19 (total n overestimate SD)
      pm3D <- rnorm(1,1.46, 0.12*sqrt(20)/1.96)    # 95% CI 0.12
      px2D <- rnorm(1,-0.79, 0.09*sqrt(123)/1.96)  # 95% CI ± 0.09            n = 123
      px3D <- rnorm(1,-0.86, 0.06*sqrt(131)/1.96)  # 95% CI ± 0.06            n = 131
      al2D <- rnorm(1,0.68, 0.12*sqrt(127)/1.96)   # CI 0.12 n = 127 (we use the total n=255 / 2 )
      al3D <- rnorm(1,1.05, 0.08*sqrt(128)/1.96)   # CI 0.08 n = 128
      d2D <- filter(da, {{int_dim}}=="2D") %>%
        mutate(mR=if_else({{res_mm}} < 0, 10^int_pm2D*{{con_mm}}^pm2D , {{res_mm}} ),
               xR=if_else({{res_den}} < 0 , 10^x02D * mR^px2D,{{res_den}}),
               alfa=10^int_al2D *{{con_mm}}^al2D,
               qRC = alfa*xR*mR/{{con_mm}})

      d3D <- filter(da, {{int_dim}}=="3D") %>%
        mutate(mR=if_else({{res_mm}} < 0, 10^int_pm3D*{{con_mm}}^pm3D , {{res_mm}} ),
               xR=if_else({{res_den}} < 0 , 10^x03D * mR^px3D,{{res_den}}),
               alfa = 10^int_al3D * {{con_mm}}^al3D,
               qRC = alfa*xR*mR/{{con_mm}}
        )
      bind_rows(d2D,d3D)

    })
    da <- bind_rows(da)
  }
  return(da)
#
}


