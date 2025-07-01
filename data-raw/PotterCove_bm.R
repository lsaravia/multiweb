## code to prepare `PotterCove_bm` dataset goes here
p <- read.delim("data-raw/Potter_iD_PC.txt", dec=",")
#
# Conversion body mass and biomass from g to kg
p <- p %>% mutate(c_bodymass=c_bodymass/1000, r_bodymass=r_bodymass/1000, r_biomass=r_biomass/1000)
#
# Assign 1 to resources whose body size is decoupled from the consumer's body size
p <- p %>% mutate(r_bodymass=if_else(is.na(r_bodymass), 1, r_bodymass)) # MAA, necromass, fresh and aged detritus
p <- p %>% mutate(r_density=if_else(r_bodymass==1, 1, r_density))
#
# Manually assign 1 in r_density and r_bodymass to sponge species
p <- p %>% mutate(r_density=if_else(resource=="Chalinidae", 1, r_density))
p <- p %>% mutate(r_bodymass=if_else(resource=="Chalinidae", 1, r_bodymass))
p <- p %>% mutate(r_density=if_else(resource=="Dendrilla antarctica", 1, r_density))
p <- p %>% mutate(r_bodymass=if_else(resource=="Dendrilla antarctica", 1, r_bodymass))
p <- p %>% mutate(r_density=if_else(resource=="Porifera", 1, r_density))
p <- p %>% mutate(r_bodymass=if_else(resource=="Porifera", 1, r_bodymass))
p <- p %>% mutate(r_density=if_else(resource=="Rosella antartica", 1, r_density))
p <- p %>% mutate(r_bodymass=if_else(resource=="Rosella antartica", 1, r_bodymass))
p <- p %>% mutate(r_density=if_else(resource=="Rossella sp.", 1, r_density))
p <- p %>% mutate(r_bodymass=if_else(resource=="Rossella sp.", 1, r_bodymass))
p <- p %>% mutate(r_density=if_else(resource=="Stylo_Myca", 1, r_density))
p <- p %>% mutate(r_bodymass=if_else(resource=="Stylo_Myca", 1, r_bodymass))
#
# When resources density information is not available assign the value -999 so the function can calculate the body mass based on allometric relationship of its consumers and the resource density according Pawar et al. (2012) original equation
p$r_density[is.na(p$r_density)] <- -999
#
# Check that there is no body mass and/or density NA data for resources and consumers
# p %>% filter(is.na(r_bodymass))
# p %>% filter(is.na(c_bodymass))
# p %>% filter(is.na(r_density))
PotterCove_bm <- p
usethis::use_data(PotterCove_bm, overwrite = TRUE)
