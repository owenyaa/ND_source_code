This folder mainly contains the following files,
Main file: Main_ program_ TMRM.m
The parameters of the algorithm can be modified and set in this document, including the order of series solution retained, the initial value of iteration, etc.

Sensitivity verification document:Verification_sensitivity_matrix.m
The function of this document is to verify that the response sensitivity respect to the parameters is calculated correctly.

Residual calculation and sensitivity calculation subroutine:cal_residual.m
The function of this document is to calculate the residual vector and sensitivity matrix of the system.

The subroutines related to the ERSA algorithm are as follows:csvd.m; l_corner.m; l_curve.m; lcfun.m; plot_lc.m; tikhonov.m.
These files are subroutines that iteratively solve the optimal solution of the objective function through the ERSA algorithm.
