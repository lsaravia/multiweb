# Version  0.6.60.0000

* Add function calc_QSS_extinction_dif_grp to calculate the difference in QSS with a group of species.  


# Version  0.6.50.0000

* Updated version of calc_interaction_intensity() with different sd for the normal distributions used
  for simulations around the mean of the coeficients estimated by Pawar 2012.  

# Version  0.6.31.0000

* Updated documentation of calc_interaction_intensity function

* Add @importFrom kSamples to calc_QSS_extinction_dif

# Version  0.6.30.0000

* Updated defaults of calc_QSS(), now the parameters defaults are `negative=-10, positive=1, selfDamping=-1` 
equivalent to a ecological transfer efficieny of 10%.
Previusly they were `negative=-10, positive=.1, selfDamping=-1`.

* Updated documentation with better explication of calc_interaction_intensity function
