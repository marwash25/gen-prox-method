# gen-prox-method

Matlab code for reproducing the experimental results of the paper [General Proximal Gradient Method: A Case for Non-Euclidean Norms](https://infoscience.epfl.ch/record/230391?ln=en), by Marwa El Halabi, Ya-Ping Hsieh, Bang Vu, Quang Nguyen, and Volkan Cevher.

## Usage:

### Sparse Linear Regression experiment

To reproduce the figures (Figure 1) for this experiment run the script `main_l1reg_plot.m` (with `s = 10, 50, 100`).

To reproduce the corresponding results call the function:

` main_l1reg(p,n,s,maxit,noise_sigma)`

with the following parameters (see paper for details):
* `p` dimension
* `n` number of measurements
* `s` sparsity
* `maxit` maximum number of iterations
* `noise_sigma` noise standard deviation

**Dependencies**: CVX [GB14].

### Latent group Lasso experiment

To reproduce the time table (Table 1) for this experiment run the script `time_proxLGL_print.m` 

To reproduce the corresponding results call the function:

`testproxLGL_run(sizeGrp)`

with `sizeGrp` the size of the groups (set to 10).

To reproduce the figures (Figure 2) for this experiment run the script `main_LGLreg_Plot.m` 

To reproduce the corresponding results call the function:

`main_LGLreg_run(p,n,s,sizeGrp,randomGrps,maxit,adapt_tol,tol)`

with the following parameters (see paper for details):
* `p` dimension
* `n` number of measurements
* `s` number of active groups
* `sizeGrp` size of groups
* `randomGrps` type of groups (`1` for random groups, `0` for interval groups)
* `maxit` maximum number of iterations
* `adapt_tol` flag for prox precision (`1` to decrease linearly with iterations, `0` otherwise)
* `tol` prox initial precision

**Dependencies**: CVX [GB14], and Gurobi [G15].

## References:
 
* [G15] Gurobi Optimization, version 6.0.4. https://www.gurobi.com/, 2015.
* [GB14] M. Grant and S. Boyd. CVX: Matlab software for disciplined convex programming, version 2.1. http://cvxr.com/cvx/, Mar. 2014.


