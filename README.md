# SU2deParaguai

Finite Volume Method for the solution of the compressible Euler equations. Based on SU2, Blazek's book and Toro's book on Riemann Solvers. Currently unfinished. 
TODO: finish upwind and central fluxes, add boundary conditions, main loop and time integration. Add viscous fluxes later? Maybe.
Currently, only computes gradients of solution (Green-Gauss and WLS), Rusanov flux and outpus to .vtu. Reads in .su2 mesh files.