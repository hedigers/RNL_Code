
Robust Nonlinear Shrinkage (R-NL)

The main function of R-NL in matlab and R and the relevant code to replicate the simulations of the paper
'R-NL: Covariance Matrix Estimation for Elliptical Distributions based on Nonlinear Shrinkage',
see: https://ieeexplore.ieee.org/abstract/document/10109124

Note: 'run_sim_tdist.m' is the main file to replicate the simulations.

To run R-NL or R-C-NL (the correlation based variant) on some nxp dimensional data matrix X, run
RNL(X,0) for the classical R-NL and RNL(X) for R-C-NL (recommended). 
