#' Calculates the interaction intensity of a food web using the metabolic theory and
#' the interaction dimensionality
#'
#' The function uses the body mass in Kg of predator/consumer and prey/resources and the dimensionality of the interaction as source data,
#' then the interaction intensity is estimated with all the coefficients from Pawar (2012) as `alfa*xR*mR/mC`, where `alpha` is the search rate `xR`
#' the resource density, `mR` the resource body mass and `mC` the consumer body mass. This value of the interaction strength quantifies
#' the effect of the predator on the prey by biomass unit of the predator. Assuming a Lotka-Volterra model is equivalent to the entry A(i,j) of the community matrix, where i is the
#' prey and j the predator.
#'
#' If the resource density is unknown (parameter `res_den`) you could set the column to a less than 0 value; and it
#' will be estimated according to the equation S18 and supplementary figures 2i & j (individuals/m2 - m3)
#'
#' If the mean mass of the resource for detritus or sediment (parameter `res_mm`) is unknown, it can be designated as
#' negative. This will result in the calculation of the resource body mass (in kilograms) using allometric formulas given
#' in the Equation S9 and Supplementary Figures 2c & d from the paper. This is only valid when the size ratios tend
#' to be optimal.
#'
#' If the Biomass of the resource is known you should use it as `res_mm` and set `res_den` to 1. This is the best choice to
#' avoid the previous allometric calculations of `res_den` and  `res_mm` when they are unknown.
#'
#' If resource size `res_mm` and resource density `res_den` are decoupled from consumer size you could assign 1 to both see
#' pag 487 **Dimensionality and trophic interaction strengths** in Pawar's paper.
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
  pm2D = 0.73 # 0.67 # 95% CI 0.17 sd = SE * sqrt(n) / 1.96, n = 39 (total n overestimate SD)
  pm3D = 0.92 # 1.46 # 95% CI 0.12
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
  det <- da %>% filter({{res_mm}}<0,{{res_den}}<0)
  if( nrow(det) > 0 ) warning(paste("Resource size res_mm and resource density res_den are both < 0 for", nrow(det), "rows\n and the estimationc will be highly uncertain"))
  det <- da %>% filter({{res_mm}}==0)
  if( nrow(det) > 0 ) warning(paste("If resource size res_mm == 0 the interaction strength will be 0 for", nrow(det), "rows"))

  det <- da %>% filter({{con_mm}}<0 )
  if( nrow(det) > 0 ) stop(paste("Consumer size con_mm cannot be < 0 for", nrow(det), "rows"))


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
      pm2D <- rnorm(1,0.73, 0.10)    # We assume that the values reported are the SD or 95% CI, in that last case we may be overestimating the true sd
      pm3D <- rnorm(1,0.92, 0.08)    # 95% CI 0.12
      px2D <- rnorm(1,-0.79, 0.09)   # 95% CI ± 0.09            n = 123
      px3D <- rnorm(1,-0.86, 0.06)   # 95% CI ± 0.06            n = 131
      al2D <- rnorm(1,0.68, 0.12)    # CI 0.12 n = 127 (we use the total n=255 / 2 )
      al3D <- rnorm(1,1.05, 0.08)    # CI 0.08 n = 128
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


