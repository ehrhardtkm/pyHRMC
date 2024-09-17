Getting Started
===

Input files
---
Running a simulation with pyHRMC requires several input files that contain calculation parameters, experimental data, and an initial structure configuration. To begin, create a new directory with the following files:

- runfile.py
- your_epdf.gr
- initial_structure.vasp
- LAMMPS input files:
    - in.lmp
    - in_init.lmp
    - in_accept.lmp
- FLARE force field model (Optional)