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
Additionally, if using a potential model trained from FLARE, follow instructions in [their documentation](https://colab.research.google.com/drive/1qgGlfu1BlXQgSrnolS4c4AYeZ-2TaX5Y#scrollTo=VYJVQ7XSWUEq) to include this in the input script.
```
pair_style    flare
pair_coeff    * * lmp.flare
```
**Atom groups:** Add the mass and element names. The numbering of the atom types must follow the order of atoms in the original structure file.
```
mass {atom_type} {float}
group {element_name} type {atom_type}
```
**Compute and outputs:** Include commands to compute and print any additional desired variables. Any commands in this section are optional. The program will automatically extract the potential energy of the structure and dump the trajectory into the XDATCAR file, but additional desired outputs can be included here and dumped into a separate file.

If using a FLARE potential, the uncertainty of the potential can be dumped on a per-atom basis. If the uncertainty is computed, LAMMPS can dump the trajectory in a file which can be used to visualize the uncertainty in OVITO. Detailed instructions can be found in this [tutorial notebook](https://colab.research.google.com/drive/1qgGlfu1BlXQgSrnolS4c4AYeZ-2TaX5Y#scrollTo=TH2-qLxsaN89). In the input script, include the following lines to compute the uncertainty-per-atom and the maximum atomic uncertainty and dump the uncertainty to a trajectory file called `output.dump`.

```
compute unc all flare/std/atom L_inv_lmp.flare sparse_desc_lmp.flare
compute MaxUnc all reduce max c_unc

# Define the dump file and append to it
dump dump_all all custom 1 output.dump id type x y z c_unc
dump_modify dump_all append yes sort id
```

While a thermo and thermo_style command can also be included, it is not very useful, given that the LAMMPS is called repeatedly in HRMC instead of running a single LAMMPS run over the entire simulation. 

**Run:** Run the LAMMPS code for 0 timesteps, meaning that it will only perform the desired computes and not perform any molecular dynamics
```
# Run the simulation
run 0
```

### Example
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