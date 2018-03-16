## Example simulation of "parasite-induced migratory stalling"

Simulation shows the density along a migration corrdior from 0 to 2000 km of 

  (a) host population *N(x,t)*, 
  
  (b) the mean number of parasites per host *m(x,t)*, 
  
  (c) the variance-to-mean ratio (VMR) of the parasite distribution *A(x,t)*, and 
  
  (d) free-living parasite larvae, *L(x,t)*. 
  
Solid lines correspond to moving hosts and parasites of moving hosts and dashed lines correspond to stationary hosts and parasites of stationary hosts (note that free-living larvae are only stationary, as it is assumed that their movement is negligible relative to the migration distance of hosts). 

![Alt Text](https://github.com/sjpeacock/Migration_model/blob/master/Stalling/StallingGIF.gif)

The animation moves through time as hosts migrate from left to right. Initially, all hosts begin moving at constant speed, but there is an "infection hotspot" around 750 km with a peak in density of free-lviing parasites in the environment. Hosts stop moving and become stationary at parasite-mediated rate \omega + \theta \times \hat{m}(x,t). As hosts hit the "infection hotspot", their parasite burden increases (b), and an increasing number of hosts stop their migration. This results in a further build up of free-living parasites, higher and higher parasite burdens of stationary hosts at that spot, and enventual "migratory stalling" of hosts. The final distribution of hosts at the "end" of the migration is highest at the infection hotspot around 750 km, and practically no hosts have made it to the end of the migration corridor at 2000 km.

Parameters for the simulation are: \beta = 0, \mu = 0.1, \alpha = 0.1, \sigma = 5, \omega = 1, \gamma = 1, \theta = 2, and c = 10000, \kappa = 10, \lambda = 0.8, \mu_L = 5.
