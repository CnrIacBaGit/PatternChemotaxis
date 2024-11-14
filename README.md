# PatternChemotaxis
Angela Monti Institute for Applied Mathematics (IAC), CNR, Bari, Italy. Mail: a.monti@ba.iac.cnr.it

Pattern_chemotaxis is a MATLAB routine (version R2024b) that produces spatial patterns induced by chemotaxis instability.

Models:

- The MOMOS model, that gives rise to stripes or spots;
- The Mimura-Tsujikawa model, that gives rise to hexagons.

The repository contains:

- the scripts "bif_diagram_ " that produce the bifurcation diagram for chemotaxis driven instability for the models;

- the scripts "spatial_discretization_" that generate the discretization of the Laplace operator and of the chemotactic term with homogeneous Neumann or periodic boundary conditions;

- Scripts for generating data (solve the models with Symplectic schemes)
  - the script "Mimura_hexagons_nld_neumann", that computes the solution and generates the plot of the solution u at the final time and the plots of the spatial mean and error between two consecutive solutions;
  - the script "MOMOS_stripes_nld_pbcs", that computes the solution and generates the plot of the solution u at the final time and the plots of the spatial mean and error between two consecutive solutions;
  - the script "MOMOS_spots_nld_pbcs", that computes the solution and generates the plot of the solution u at the final time and the plots of the spatial mean and error between two consecutive solutions;
  - the function "pDMD" (by Alla et al, 2024) for the application of the piecewise DMD, that recalls "rqb", "rdmd_qb" and "new_reconstruction"
  - the scripts "pDMD_ " that recall the scripts for generating data and apply the pDMD, to generate the plot of the reconstructed solution at the final time and the plot of the error over time.

- the script "future_prediction" for the application of the pDMD algorithm to reconstruct data in [0,tbar] and make future prediction in (tbar,T]. This script recalls the function "fun_pdmd_prediction"
- the folder "Upwind_codes" for generating data by using second order upwind schemes for the discretization of the chemotaxis

The routine has been implemented and developed by Angela Monti. It can be used under the conditions of CC-BY-NC 2.0

A full description of the model is available in: A. Monti, F. Diele, D. Lacitignola, C. Marangi, Patterns in soil organic carbon dynamics: integrating microbial activity, chemotaxis, and data-driven approaches (submitted), 2024
