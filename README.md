# PatternChemotaxis

Angela Monti
Institute for Applied Mathematics (IAC), CNR, Bari, Italy. 
Mail: a.monti@ba.iac.cnr.it

Pattern_chemotaxis is a MATLAB routine (version R2023a) that produces spatial patterns induced by chemotaxis instability. 

Models:
- The MOMOS model, that gives rise to stripes or spots;
- The Mimura-Tsujikawa model, that gives rise to hexagons.

The repository contains:
- the scripts "bif_diagram_ " that produce the bifurcation diagram for chemotaxis driven instability for the models;

- the script "spatial_discretization" that generates the discretization of the Laplace operator and of the chemotactic term;

- Scripts for generating data (solve the models with Symplectic schemes)
1. the script Mimura_hexagons, that computes the solution and generates the plot of the solution u at the final time and the plot of the spatial mean;

2. the script MOMOS_stripes, that computes the solution and generates the plot of the solution u at the final time and the plot of the spatial mean;

3. the script MOMOS_spots, that computes the solution and generates the plot of the solution u at the final time and the plot of the spatial mean;

- the function pDMD (by Alla et al, 2024) for the application of the piecewise DMD, that recalls rqb, rdmd_qb and new_reconstruction

- the scripts "pDMD_  " that recall the scripts for generating data and apply the pDMD, to generate the plot of the reconstructed solution at the final time and the plot of the error over time. 

The routine has been implemented and developed by Angela Monti. It can be used under the conditions of CC-BY-NC 2.0

A full description of the model is available in: 
A. Monti, F. Diele, D. Lacitignola, C. Marangi, Patterns in soil organic carbon dynamics: integrating microbial activity, chemotaxis, and data-driven approaches (submitted), 2024

