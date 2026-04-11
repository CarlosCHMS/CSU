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

The solver considers values in the international system of units.

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

Several commands can be understood from the case folders provided. Here, one of the `input.dat` files is described in detail:

- **pressure, 1e5**: pressure at the inlet

- **mach, 0.1**: Mach number at the inlet

- **temperature, 300.**: temperature at the inlet

- **nx, 1.**: x-component of the unit velocity vector at the inlet 

- **ny, 0.**: y-component of the unit velocity vector at the inlet

- **pout, 1e5**: pressure at the outlet

- **BC:wall, wall**: wall boundary associated with the wall boundary condition

- **BC:outlet, outlet**: outlet boundary associated with the outlet boundary condition

- **BC:sym, symmetry**: symmetry boundary associated with the symmetry boundary condition

- **BC:inlet, inlet**: inlet boundary associated with the inlet boundary condition

- **BC:boundary_name, boundary_type**: general format for defining boundary conditions

- **Nmax, 40000**: maximum number of iterations

- **order, 2**: order of the solver

- **threads, 4**: number of threads for parallel processing

- **flux, AUSMpup2**: flux scheme type. Available options are: ROE, AUSM, AUSMDV, AUSMpup, and AUSMpup2

- **axisymmetric, 0**: axisymmetric model flag disabled

- **CFL, 1e3**: Courant number. Must be smaller than 1.0 for explicit solvers (RK)

- **laminar, 0**: laminar model flag enabled

- **sa, 1**: Spalart–Allmaras turbulence model flag enabled (RANS)

- **sst, 0**: k-ω SST turbulence model flag disabled (RANS)

- **restart, 0**: restart from a previous solution disabled

- **limK, 0.1**: K factor in the Venkatakrishnan limiter

- **wImp, 1.5**: relaxation factor in the LUSGS implicit scheme

- **timeScheme, LUSGS**: time integration scheme. Options are: RK (Runge–Kutta), LUSGS

- **TP, 0**: thermally perfect gas model flag disabled

## Boundary condition types

The following boundary types are available:

- **wall**: adiabatic wall
- **wallT**: isothermal wall
- **inlet**: subsonic or supersonic inlet
- **outlet**: subsonic or supersonic outlet
- **symmetry**: symmetry boundary

For the **wallT** condition, the command:

- **Twall, "temperature"**

must be included in `input.dat` to specify the wall temperature ("temperature" can be any numeric value).

