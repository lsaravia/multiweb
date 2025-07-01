# Version  0.8.50.0000

* New function plot_troph_level_ggraph and new dataset PotterCove_bm


# Version  0.8.40.0000

* Updated calc_QSS to effectively use the parameter selfDamping, before this changes it was used only for positive values of 
the diagonal of the adyacency matrix. Now the equivalent to the previous version is to set selfDamping=0.

* Updated toGLVadjMat which does not take into account the names of the networks in the mgraph object, now it uses the position
  the negative/positive/trophic.

# Version  0.8.30.0000

* add aggregate_multiplex_network() function to aggregate a multiplex network into a single layer network

# Version  0.8.20.0000

* add plot_multiplex_modules() function to plot the modules of a multiplex network
 
# Version  0.8.10.0000

* Improve plot_multi3D() and added harmonized_node_sets to add all the nodes to all layers, needed for plot_multi3D().

# Version  0.8.00.0000

* Added functions

  - plot_multi3D() was adapted from package muxViz.
  - convert_infomap_result_to_muxviz() to convert the output of `run_infomap_multi()` to 
    muxViz-Compatible Format.

# Version  0.7.90.0000

* The function plot_troph_level_ggplot() was updated to use return the ggplot2 object.

# Version  0.7.80.0000

* Added functions 
  - run_infomap_multi to run the python implementation of Infomap for node aligned multiplex networks
  - Other aux functions like write_multilayer_network() convert_to_intra_format()

# Version  0.7.40.0000

* Added functions 
  - run_infomap to run the python implementation of Infomap 
  - calc_centrality instead of calc_eigencentrality to calculate the centrality of a network using any defined measure


# Version  0.7.30.0000

* Added functions calc_svd_entropy shuffle_network_deg shuffle_network_deg_svd calc_eigencentrality calc_svd_entropy_importance

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
