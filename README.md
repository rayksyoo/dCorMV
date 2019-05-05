# dCorMV
Calculate multivariate distance correlation (Reference: Yoo et al. Multivariate approaches improve the reliability and validity of functional connectivity and prediction of individual behaviors. NeuroImage. 2019 )

This is to calculate multivariate distance correlation.
Input X is 1D cell array of length m. Each cell contains 2D (n by p) time-series matrices with n time points and p voxels. Here, p can be different across m, but n should be the same.

The main output is m by m distance correlation matrix.
