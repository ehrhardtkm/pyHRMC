LAMMPS input and potential files
===

In pyHRMC, the structural energy is computed for every atom hop with LAMMPS. LAMMPS requires input files that act as a "to-do list" for the program, including information about the structure, desired computes, and running LAMMPS. For more detailed documentation, please refer to the official [LAMMPS documentation](https://docs.lammps.org/Manual.html). 

Curently, pyHRMC requires three input files. In the future, this wil be modified to reduce the number of files to one. 

in.lmp
---
Four sections are required for the in.lmp file: Initialization, Atom groups, Compute and outputs, and Run. 

**Initialization:** This is used to select the units for LAMMPS to use, read in the structure, and specify potential parameters. Because pyHRMC will automatically write a structure data file in the LAMMPS format, always include the following line:
```
read_data data.lmp
```
Additionally, is using a potential model trained from FLARE, follow instructions in [their documentation](https://colab.research.google.com/drive/1qgGlfu1BlXQgSrnolS4c4AYeZ-2TaX5Y#scrollTo=VYJVQ7XSWUEq) to include this in the input script.
```
pair_style    flare
pair_coeff    * * lmp.flare
```
**Atom groups:** Add the mass and element names of the elements in the cell. The numbering of the atom types will follow the order of atoms in the original structure file
```
mass {atom_type} {float}

group {element_name} type {atom_type}
```
**Compute and outputs**
```
thermo 1
thermo_style custom step pe ke etotal
```

**Run:** Run the LAMMPS code for 0 timesteps, meaning that it will only perform the desired computes and not perform any molecular dynamics
```
# Run the simulation
run 0
write_data out4.txt nocoeff nofix
```

Example
---
```
# Initialization
atom_style atomic
units metal
boundary p p p
atom_modify sort 0 0.0
newton on

read_data data.lmp

### interactions
pair_style flare
pair_coeff * * s_lmp.flare

# Atom groups
mass 1 26.981539
mass 2 15.9994

group Al type 1
group O type 2


#Compute and outputs
compute unc all flare/std/atom s_L_inv_lmp.flare s_sparse_desc_lmp.flare
compute MaxUnc all reduce max c_unc

thermo 1
thermo_style custom step pe ke etotal


# Run the simulation
run 0
write_data out4.txt nocoeff nofix
```

in_init.lmp
---
This script requires additional lines to save each step in a file called "last_step.txt" to allow LAMMPS to restart at the appropriate timestep. Add the following lines before the run command in the original in.lmp script:
```
# Save the final timestep to a file
variable    last_step equal "step"
shell       echo ${last_step} > last_step.txt
```

in_accept.lmp
---
This script requires an additional line to read in the current timestep of the simulation. Add the following line before the run command in the original in.lmp script:

```
# Read the last timestep from the file
reset_timestep ${last_step}
```