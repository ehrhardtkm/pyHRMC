Getting Started
===

## List of input files

Running a simulation with pyHRMC requires several input files that contain calculation parameters, experimental data, and an initial structure configuration. To begin, create a new directory with the following files:

- runfile.py
- your_epdf.gr
- initial_structure.vasp
- LAMMPS input files:
    - in.lmp
    - in_init.lmp
    - in_accept.lmp
- FLARE force field model (optional)


## your_epdf.gr

The experimental ePDF data to be used in the simulation must be prepared as a file with tabulated r and G(r) values that are separated by a comma. For example, this file might look like:

```
r,gr
0.0125, 0
0.025, 0
0.0375, 0
0.05, 0
0.0625, 0
...
...
...
9.9625, -0.0011
9.975, -0.00183
9.9875, -0.00237
10, -0.00271
```

## initial_structure.vasp

A structure file must be prepared in the POSCAR/VASP format for the HRMC simulation to use at a starting configuration of atoms. For a slab simulation, sufficient vacuum space must be inserted above and below the slab to simulate these surfaces. 

The initial structure file must include the correct number of atoms and stoichiometry that are to be simulated, since HRMC does not affect the composition of the simulated cell. Additionally, if performing a bulk simulation, the cell must also be configured with the correct density. In slab simulations, density may decrease to a minor degree, but the slab will not expand to fill the vacuum space in the z direction.