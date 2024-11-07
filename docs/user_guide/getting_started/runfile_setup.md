Runfile setup
===
All calculation parameters are defined in runfile.py, which instantiates the RMC class before calling the run_rmc() method to perform the simulation. A runfile should begin with necessary imports:
```
from pyHRMC.core.rmc import RMC
import multiprocessing
```
## Creating RMC object

To run a simulation, the user must first create a RMC class object with the required arguments. These arguments determine the quench rate of the error and energy constraints.
```
rmc = RMC(
    experimental_G_csv=str(), 
    sigma=float(), 
    q_scatter=float(), 
    q_temp=float(), 
    init_temp=float(), 
    dump_freq = 5000
    )
```
- **experimental_G_csv:**  filename of experiment ePDF
- **sigma:** initial value of sigma, the value is used to calculate the probability of accepting a step that increases the error (i.e. is the worse fit to experimental ePDF). To allow proper convergence, this value should **not** start at zero, which means that steps that increase error are never accepted. Using `sigma=0` for the entire simulation will prevent the structure sufficiently exploring the structural space to find a solution with a low error. For details, see [Background and Theory](https://ehrhardtkm.github.io/pyHRMC/background_theory/)
- **q_scatter:** quench rate for error constraint (rate of decreasing sigma). If a value of q_scatter < 1 is used, `sigma` will gradually decrease over the course of the simulation.
- **q_temp:**  quench rate for energy constraint (rate of decreasing temperature). If a value of q_temp < 1 is used, this value will gradually decrease over the course of the simulation.
- **init_temp:** initial temperature, which is used to calculate the Boltzmann-type probability of accepting an energetically unfavorable step. The program will not quench the temperature below 100 K, so an `init_temp` below this value will remain unchanged regardless of the `q_temp`. For details, see [Background and Theory](https://ehrhardtkm.github.io/pyHRMC/background_theory/)
- **dump_freq:** the frequency at which to update the [pdf.txt and pdfs.png files](https://ehrhardtkm.github.io/pyHRMC/user_guide/getting_started/outputs), as well as output the [maximum atomic uncertainty](https://ehrhardtkm.github.io/pyHRMC/user_guide/getting_started/lammps_input) as computed by LAMMPS from a FLARE potential. Defaults to every 5000 HRMC iterations

The choice of `sigma`, `init_temp`, and their corresponding quenching rates `q_scatter` and `q_temp` is highly important to allow the simulation to properly converge. Many examples of simulations and their parameters used can be found in the literature, some of which have been listed in the [references section](https://ehrhardtkm.github.io/pyHRMC/background_theory/).

## Calling the run_rmc() method

To perform a simulation, use the run_rmc() method and further specify simulation parameters and constraints.
```
rmc.run_rmc(
    num_processes=int(), 
    initial_structure=str(), 
    keV=int(), 
    prdf_args={“bin_size:float}, 
    transformations={“AtomHop”:{}}, 
    validators=dict(),
    charges=dict(),
    pdf_cutoff=float(),
    max_steps=int()
    )
```
- **num_processes:** the number of processes to use in the simulation. If a serial simulation is desired, set to 1. If greater than 1, pyHRMC performs the simulation in iterations. In each iteration, a number of Monte Carlo steps equal to `num_processes` will be independently generated and validated in parallel. Accepted individual steps will then be merged and re-tested in a single structure, which constitutes one iteration. Therefore, the total number of steps will always be equal to `num_processes` x the number of completed iterations. This distinction is particularly important when selecting `q_scatter` and `q_temp` and interpreting output files.
- **initial_structure:** the filename of the structure file to use as an initial configuration. This file must be in a POSCAR/VASP format. The structure should have the correct stoichiometry, sincve HRMC does not affect the composition of the structure. 
- **keV:** the accelerating voltage of the electron microsope used to capture scattering data
- **prdf_args**
    - **bin_size:** the bin size of the computed PDFs, in angstroms. Default is 0.04 angstroms.
- **transformations**
    - **AtomHop:** currently the only kind of transformation available in pyHRMC. This transformation creates a spherical grid of radius 1 Angstrom around an atom and selects a cell within the grid randomly. The coordinates of the selected cell are treated as the new coordinates of the atom. This transformer only moves one atom at a time, but this action can be parallelized for major speedups and yield structures with multiple moved atoms in a single step.
- **validators:** details on available validators can be found on the [validator selection](http://127.0.0.1:8000/user_guide/validator_selection/) page.
    - **SlabThickness**
    - **Coordination** 
        - **BulkCoordinationRange**
        - **SurfaceCoordinationRange**
        - **SurfaceDistance**
    - **DistancesCoordination**
        - **MinDistances**
        - **BulkCoordinationRange**
        - **SurfaceCoordinationRange**
        - **SurfaceDistance**
    - **SiteDistances**
    - **TargetDensity**
- **charges:** this is an optional argument that the user may use to bypass the built-in method of interpolating partial charges for each of the elements in the structure. Additionally, if users desire that no partial charges are used and neutral atoms are assumed, they may set the charge values to 0. If using this argument, follow the following format:
'''
charges = {
    "species1": float(charge),
    "species2" = float(charge),
    ...
    }
'''
- **pdf_cutoff:** this is the distance, in Angstroms, below which will not be included in calculating the PDF error. The value will default to 1.6 Angstroms unless otherwise specified by the user.
- **gaussian_blur:** this parameter is used when calculating the structure's PDF to determine the degree of Guassian smearing that is used from [scipy.ndimage.gaussian_filter1d](https://docs.scipy.org/doc/scipy/reference/generated/scipy.ndimage.gaussian_filter1d.html). This value corresponds to sigma in the linked scipy documentation, but is included here as `gaussian_blur` to avoid confusion.
- **max_steps:** maximum number of steps to perform in the simulation. If running in parallel, note that all atom transformations performed in parallel will be considered as a single step. 

## Example
An example runfile.py for a simulation of amorphous alumina is shown below:
```
from pyhrmc.core.rmc import RMC

if __name__ == "__main__":
   
    rmc = RMC(experimental_G_csv="gr.txt", sigma = 0.2, q_scatter = 1, q_temp = 0.99995, init_temp = 1000, dump_freq = 500)
    rmc.run_rmc(
        num_processes = 16,
        initial_structure="structure.vasp",
        experimental_G_csv="al2o3_5nm_gr.txt",
        keV=200,
        prdf_args={"bin_size": 0.04},
        # error_constant=0.05,
        transformations={
            # "ASECoordinatePerturbation":{}, #simmate version of RattleMutation from ase, rattle_prob auto set to 0.3
            "AtomHop": {}, # consider adding a second, smaller step size "max_step": 0.2
        },

        validators={
            "SlabThickness": {"max_thickness": 51.5},
            "Coordination": {
                "BulkCoordinationRange": {"Al": {"Al" : [0, 1], "O": [2, 7]}, "O": {"Al": [1, 4], "O": [0, 1]} },
                "SurfaceCoordinationRange": {"Al": {"Al": [0, 0], "O": [3, 7]}, "O": {"Al": [1, 4], "O": [0, 1]} },
                "SurfaceDistance": 3
                }
            },
            gaussian_blur = 2,
        max_steps= 1000000,
    )
```