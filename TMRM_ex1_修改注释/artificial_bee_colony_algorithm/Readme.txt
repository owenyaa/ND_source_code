This folder mainly contains the following files,
Main file: Main_ program_ TMRM.m
The parameters of the algorithm can be modified and set in this document, including the order of series solution retained, the initial value of iteration, etc.

Fitness function document:fitness1.m
Fitness function, calculate fitness value.

Random function: suiji.m
through random number to generate new individuals

The subroutines related to the ERSA algorithm are as follows:plot_the_response.m

For comparison, we also solved problem \eqref{eq4.10} in Case IV by the ABC algorithm. ABC algorithm found the optimal solution through 2772 iterations (the maximum number of cycles is 6000), but the corresponding to $\Vert \bm{R'}\Vert^2=25.19$. 

The final result is saved in matlab.mat .
