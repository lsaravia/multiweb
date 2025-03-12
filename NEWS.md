# Version  0.7.10.0000

* Added function generate_niche to generate a niche model with a given number of species and connectance


# Version  0.7.01.0000

* Added function plot_troph_level_ggplot() with improved network visualization


# Version  0.7.00.0000

* Added function calc_interaction_intensity2() with estimation of interaction strength using O'Gorman 2010 method


# Version  0.6.9.0000

* Updated function calc_interaction_intensity() with more error checking and warnings


# Version  0.6.8.0000

* Updated function calc_QSS_extinctions_seq() and calc_QSS_extinction_dif_grp() to return the raw values of the QSS simulations, KA and AD test are not performed.

# Version  0.6.7.0000

* Updated function calc_QSS_extinction_dif() to return the raw values of the QSS simulations and the difference between them, KA and AD test are not performed, the null model with mean IS is also deleted.


# Version  0.6.67.0000

* Updated function calc_topological_roles() to accept a community object

# Version  0.6.66.0000

* Updated function classify_topological_roles() when all species are module specialists 'modspe'


# Version  0.6.65.0000

* Updated README with function fromIgraphToMgraph

* Fixed calc_QSS for mgraph objects 

# Version  0.6.61.0000

* Fixed warmings of function calc_topological_roles. 

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
