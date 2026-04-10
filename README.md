# CSU

A CFD code for unstructured meshs.

To run the code use the command (linux): $ python3 CSU.py ./caseFolder/

This command will compile the c code, run the case and display the results.

# Theory

The code is a density-based finite volume CFD solver.

More details about the algorithm are presented in the paper: https://jatm.com.br/jatm/article/view/1317/989

Since the publication of this paper, several improvements have been made, including:

-The k-omega SST turbulence model;

-Flux algorithms: AUSM, AUSM+up, and AUSM+up2;

-The implicit LUSGS algorithm.
