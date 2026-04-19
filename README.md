# CSU

A CFD code for unstructured meshs.

To run the code use the command (linux): 

```bash
python CSU.py ./caseFolder/
```

This command will compile the c code, run the case and display the results.

# Theory

The code is a density-based finite volume CFD solver. It supports 2-D planar and axisymmetric meshes.

More details about the algorithm are presented in the paper: https://jatm.com.br/jatm/article/view/1317/989

Since the publication of this paper, several improvements have been made, including:

- The k-omega SST turbulence model;

- Flux algorithms: AUSM, AUSM+, AUSM+up, and AUSM+up2;

- The implicit LUSGS algorithm;

- The thermally perfect gas model.

The solver considers values in the international system of units.

# How it works

The cases are organized into folders. Inside each case folder, there are:

- **input.ini**: contains the commands that will be used to run CSU;
- **mesh.su2**: the mesh file. It uses the same format as the SU2 software and can be generated using Gmsh;
- **analysis.py**: a Python script that is executed after the simulation to display the results.

There are also:

- **geometry.geo**: an input file for Gmsh used to generate the mesh;

All solution data will be saved inside the case folder. You can run the `analysis.py` script from within the case folder using the following command:

```bash
python ./analysis.py ./
```

With this command, you can view the results again without rerunning the case.

# Commands of input.ini file

The case folders provided include `input.ini` files with several types of commands that can be used in CSU. These `input.ini` files also contain comments for each command, helping to better understand how each one is used.

## Boundary condition types

The following boundary types are available:

- **wall**: adiabatic wall
- **wallT**: isothermal wall
- **inlet**: subsonic or supersonic inlet
- **outlet**: subsonic or supersonic outlet
- **symmetry**: symmetry boundary

For the **wallT** condition, the command:

- **Twall, "temperature"**

must be included in `input.ini` to specify the wall temperature ("temperature" can be any numeric value).

