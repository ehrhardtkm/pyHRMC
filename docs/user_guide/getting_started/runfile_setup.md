Runfile setup
---
All calculation parameters are defined in runfile.py, which instantiates the RMC class before calling the run_rmc() method to perform the simulation. A runfile should begin with necessary imports:
```
from pyHRMC.core.rmc import RMC
import multiprocessing
```
To run a simulation, the user must first instantiate the RMC class with the required arguments. These arguments determine the quench rate of the error and energy constraints.
```
rmc = RMC(experimental_G_csv=str, sigma=float, q_scatter=float, q_temp=float, init_temp=float)
```
- experimental_G_csv = filename of experiment ePDF
- sigma = initial value of sigma, the measure of uncertainty in ePDF data
- q_scatter = quench rate for error constraint (rate of decreasing sigma)
- q_temp =  quench rate for energy constraint (rate of decreasing temperature)
- init_temp = initial temperature, which determines the acceptance probability based on energy constraint

To perform a simulation, use the run_rmc() method and further specify simulation parameters and constraints.
```
rmc.run_rmc(num_processes=int, initial_structure=str, keV=int, prdf_args={“bin_size:float}, transformations={“AtomHop”:{}}, validators={args}, max_steps=int)
```
- num_processes = the number of processes to use in the simulation. If a serial simulation is desired, set to 1. 
- initial_structure = the filename of the structure file to use as an initial configuration. This file must be in a POSCAR/VASP format. The structure should have the correct stoichiometry, sincve HRMC does not affect the composition of the structure. 
- keV = the accelerating voltage of the electron microsope used to capture scattering data
- prdf_args
    - bin_size: the bin size of the computed PDFs, in angstroms. Default is 0.04 angstroms.
- transformations
    - AtomHop: currently the only kind of transformation available in pyHRMC. This transformation creates a spherical grid of radius 1 Angstrom around an atom and selects a cell within the grid randomly. The coordinates of the selected cell are treated as the new coordinates of the atom. This transformer only moves one atom at a time, but this action can be parallelized for major speedups and yield structures with multiple moved atoms in a single step.
- validators = details on available validators can be found on the [Validator selection](http://127.0.0.1:8000/user_guide/validator_selection/) page.
    - SlabThickness
    - Coordination 
        - BulkCoordinationRange
        - SurfaceCoordinationRange
        - SurfaceDistance

- max_steps = maximum number of steps to perform in the simulation. If running in parallel, note that all atom transformations performed in parallel will be considered as a single step. 

Example
---
An example runfile.py for a simulation of amorphous Al2O3 is shown below:
```
from pyHRMC.core.rmc import RMC
import multiprocessing

if __name__ == "__main__":

    rmc = RMC(experimental_G_csv="al2o3_5nm_gr.txt", sigma = 0.05, q_scatter = 1, q_temp = 0.999999995, init_temp = 1000)
    rmc.run_rmc(
        num_processes = 8,
        initial_structure="POSCAR_47A-d33",
        keV=200,
        prdf_args={"bin_size": 0.04},
        transformations={
            "AtomHop": {}, # consider adding a second, smaller step size "max_step": 0.2
        },

        validators={
            "SlabThickness": {"max_thickness": 47.4},
            "Coordination": {
                "BulkCoordinationRange": {"Al": {"Al" : [0, 1], "O": [2, 7]}, "O": {"Al": [1, 4], "O": [0, 1]} },
                "SurfaceCoordinationRange": {"Al": {"Al": [0, 1], "O": [2, 7]}, "O": {"Al": [1, 4], "O": [0, 1]} },
                "SurfaceDistance": 5
                }
            },

        max_steps= 1000000,
    )
```