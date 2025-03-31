Runfile setup
===
All calculation parameters are defined in runfile.py, which instantiates the RMC class before calling the run_rmc() method to perform the simulation. A runfile should begin with necessary imports:
```
from pyHRMC.core.rmc import RMC
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
    qmin=float(),
    charges=dict(),
    TCS=dict()
    pdf_cutoff=float(),
    gaussian_blur=int(),
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
- **qmin:** minimum q value of the experimental scattering data.
- **charges:** NOTE: At this time, we do not advise users to use this argument for their simulations.
- **TCS:** if desired, users my choose to use their own total scattering cross-section values. This allows HRMC calculations with X-ray or neutron scattering cross-sections, if desired. By choosing to omit this argument, the calculation will default to electron scattering cross-sections. When entering TCS values, use the following format:
```
TCS = {
    'Al': 0.011,
    'O': 0.0045
}
```
- **pdf_cutoff:** this is the distance, in Angstroms, below which will not be included in calculating the PDF error. The value will default to 1.6 Angstroms unless otherwise specified by the user.
- **gaussian_blur:** this parameter is used when calculating the structure's PDF to determine the degree of Guassian smearing that is used from [scipy.ndimage.gaussian_filter1d](https://docs.scipy.org/doc/scipy/reference/generated/scipy.ndimage.gaussian_filter1d.html). This value corresponds to sigma in the linked scipy documentation, but is included here as `gaussian_blur` to avoid confusion.
- **max_steps:** maximum number of steps to perform in the simulation. If running in parallel, note that all atom transformations performed in parallel will be considered as a single step. 

## Example
An example runfile.py for a simulation of amorphous alumina is shown below:
```
from pyhrmc.core.rmc import RMC

if __name__ == "__main__":

    rmc = RMC(experimental_G_csv="7ps_blur2.txt", sigma = 0.5, q_scatter = 0.9995, q_temp = 0.9995, init_temp = 2000, dump_freq = 500)
    rmc.run_rmc(
        num_processes = 16,
        initial_structure="crystal_d3_slab.vasp",
        experimental_G_csv="7ps_blur2.txt",
        keV=200,
        prdf_args={"bin_size": 0.04},
        transformations={
            "AtomHop": {}
        },

        validators={
            "SlabThickness": {"max_thickness": 51.5},
            "DistancesCoordination": {
                "MinDistances": {
                    ("Al", "Al"): 2.0,
                    ("Al", "O"): 1.6,
                    ("O", "O"): 2.0
                    },
                "BulkCoordinationRange": {"Al": {"Al" : [0, 0], "O": [4, 5]}, "O": {"Al": [3, 4], "O": [0, 0]} },
                "SurfaceCoordinationRange": {"Al": {"Al": [0, 0], "O": [4, 5]}, "O": {"Al": [2, 4], "O": [0, 0]} },
                "SurfaceDistance": 3
                    }
            },
        qmin = 2,
        gaussian_blur = 2,
        charges={
            "Al": 0,
            "O": 0,
            },
        max_steps= 1000000,
    )

```