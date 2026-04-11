# CSU

A CFD code for unstructured meshs.

To run the code use the command (linux): 

```bash
python3 CSU.py ./caseFolder/
```

This command will compile the c code, run the case and display the results.

# Theory

The code is a density-based finite volume CFD solver.

More details about the algorithm are presented in the paper: https://jatm.com.br/jatm/article/view/1317/989

Since the publication of this paper, several improvements have been made, including:

- The k-omega SST turbulence model;

- Flux algorithms: AUSM, AUSM+up, and AUSM+up2;

- The implicit LUSGS algorithm;

- The thermally perfect gas model.

# How it works

The cases are organized into folders. Inside each case folder, there are:

- **input.dat**: contains the commands that will be used to run CSU;
- **mesh.su2**: the mesh file. It uses the same format as the SU2 software and can be generated using Gmsh;
- **analysis.py**: a Python script that is executed after the simulation to display the results.

There are also:

- **geometry.geo**: an input file for Gmsh used to generate the mesh;
- **su2MeshReader.py**: a module used by `analysis.py` to read the `mesh.su2` file.

All solution data will be saved inside the case folder. You can run the `analysis.py` script from within the case folder using the following command:

```bash
python3 ./analysis.py ./
```

With this command, you can view the results again without rerunning the case.

# Commands of input.dat file

Several commands can be obtained from case folders presented. However, some of than must be detailed:

The types of boundaries are associated to the mesh boundaries using the commands:

BC:boundary_name, boundary_type

There are the boundary types: 

- wall: wall adiabatic
- wallT: wall isothermic
- inlet: inlet subsonic or supersonic
- outlet: outlet subsonic or supersonic
- symmetry: symmetry

