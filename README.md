# Discrete and Continuous Difference of Submodular Minimization

Code to reproduce results of the paper [Discrete and Continuous Difference of Submodular Minimization](https://arxiv.org/pdf/2506.07952).

## To reproduce results in the paper

### Integer least squares 
- run script Least-Squares/ilsq.sh for experiments with varying measurements, or Least-Squares/noisy_ilsq.sh for experiments with varying noise.
- plot results using Least-Squares/plot_measurements.m for experiments with varying measurements, or Least-Squares/plot_noisy.m for experiments with varying noise.

### Integer compressed sensing 
- run  script CS/integer_cs.sh for experiments with varying measurements, or CS/noisy_integer_cs.sh for experiments with varying noise.
- plot results using CS/plot_measurements.m for experiments with varying measurements, or CS/plot_noisy.m for experiments with varying noise.

# Citation
```
@InProceedings{orfanides2025discrete,
      title={Discrete and Continuous Difference of Submodular Minimization},
      author={George Orfanides and Tim Hoheisel and Marwa El Halabi},
      booktitle = {Proceedings of the 42nd International Conference on Machine Learning},
      year={2025},
}
```

# Acknowledgements
- We use a number of functions implemented from [1] (such as the pairwise FW algorithm for submodular minimization).  
- We use the FISTA implementation from [2].
- We use the ADMM heuristic for mixed-integer quadratic programming from [3].
- Our implementation of OMP is based on the implementation in [4].
- We use for plotting the distinguishable_colors function from [5].

[1]  Bach, Francis. "Submodular functions: from discrete to continuous domains." Mathematical Programming 175 (2019): 419-459. (http://www.di.ens.fr/~fbach/submodular_multi_online.zip) Retrieved January 22, 2024.

[2]  Beck, Amir, and Nili Guttmann-Beck. "FOM–a MATLAB toolbox of first-order methods for solving convex optimization problems." Optimization Methods and Software 34.1 (2019): 172-193. (https://www.tau.ac.il/~becka/solvers/fista) Retrieved June 4, 2024.

[3] Takapoui, Reza, et al. "A simple effective heuristic for embedded mixed-integer quadratic programming." International journal of control 93.1 (2020): 2-12. (https://github.com/cvxgrp/miqp_admm) Retrieved May 17, 2021. Retrieved December 11, 2024.

[4] Stephen Becker (2025). CoSaMP and OMP for sparse recovery (https://www.mathworks.com/matlabcentral/fileexchange/32402-cosamp-and-omp-for-sparse-recovery), MATLAB Central File Exchange. Retrieved July 6, 2025.

[5]  Holy, Tim. "Generate maximally perceptually-distinct colors." MATLAB Central File Exchange. (https://www.mathworks.com/matlabcentral/fileexchange/29702-generate-maximally-perceptually-distinct-colors) Retrieved May 17, 2021.
