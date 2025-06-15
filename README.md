# Discrete and Continuous Difference of Submodular Minimization

Code to reproduce results of the paper [Discrete and Continuous Difference of Submodular Minimization](https://arxiv.org/pdf/2506.07952)

## To reproduce results in the paper

### Integer least squares 
- run  script Least-Squares/ilsq.sh for exmperiments with varying measurements, or Least-Squares/noisy_ilsq.sh for exmperiments with varying noise.
- plot results using Least-Squares/plot_measurements.m for exmperiments with varying measurements, or Least-Squares/plot_noisy.m for exmperiments with varying noise.

### Integer compressed sensing 
- run  script CS/integer_cs.sh for exmperiments with varying measurements, or CS/noisy_integer_cs.sh for exmperiments with varying noise.
- plot results using CS/plot_measurements.m for exmperiments with varying measurements, or CS/plot_noisy.m for exmperiments with varying noise.

# Acknowledgements
- We use a number of functions implemented in [1] (such as the pairwise FW algorithm for submodular minimization). The code is available at http://www.di.ens.fr/~fbach/submodular_multi_online.zip.  
- We use the FISTA implementation from [2].
- We use for plotting the distinguishable_colors function from [3].

[1]  Bach, Francis. "Submodular functions: from discrete to continuous domains." Mathematical Programming 175 (2019): 419-459. 

[2]  Beck, Amir, and Nili Guttmann-Beck. "FOM–a MATLAB toolbox of first-order methods for solving convex optimization problems." Optimization Methods and Software 34.1 (2019): 172-193.

[3]  Holy, Tim. "Generate maximally perceptually-distinct colors." MATLAB Central File Exchange. (https://www.mathworks.com/matlabcentral/fileexchange/29702-generate-maximally-perceptually-distinct-colors) Retrieved May 17, 2021.
