User Guide
---
Running a simulation with pyHRMC requires several input files that contain calculation parameters, experimental data, and an initial structure configuration. To begin, create a new directory with the following files:
-	Runfile.py
-	Your epdf.gr
-	initial_structure.vasp
-	LAMMPS input files:
o	In.lmp
o	In_init.lmp
o	In_accept.lmp
-	FLARE force field model (Optional)
All calculation parameters are defined in runfile.py, which instantiates the RMC class before calling the run_rmc() method to perform the simulation. 
```
rmc = RMC(experimental_G_csv=str, sigma=float, q_scatter=float, q_temp=float, init_temp=float)
```
The RMC class takes arguments that determine the quench rate of the error and energy constraints.
```
Experimental_G_csv: filename of experiment ePDF
sigma: # initial value of sigma, the measure of uncertainty in ePDF data
q_scatter: # quench rate for error constraint (rate of decreasing sigma)
q_temp:  quench rate for energy constraint (rate of decreasing temperature)
init_temp: initial temperature, which determines the acceptance probability based on energy constraint
```
To perform a simulation, use the run_rmc() method and further specify simulation parameters and constraints.
```
rmc.run_rmc(num_processes=int, initial_structure=str, keV=int, prdf_args={“bin_size:float}, transformations={“AtomHop”:{}}, validators={args}, max_steps=int)
```



An example for a simulation of amorphous Al2O3 is shown below:
```
from simmate.apps.epdf_rmc.rmc import RMC
import multiprocessing

if __name__ == "__main__":

    rmc = RMC(experimental_G_csv="al2o3_5nm_gr.txt", sigma = 0.05, q_scatter = 1, q_temp = 0.999999995, init_temp = 1000)
    rmc.run_rmc(
        num_processes = 8,
        initial_structure="POSCAR_47A-d33",
        keV=200,
        prdf_args={"bin_size": 0.04},
        transformations={
            # "ASECoordinatePerturbation":{}, #simmate version of RattleMutation from ase, rattle_prob auto set to 0.3
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

Traditional RMC simulations can also be performed with pyHRMC. It’s important to note that traditional RMC is typically less accurate than HRMC, and for this reason, we don’t recommend it. However, it can be useful in some cases, such as when a suitable force field is not available for energy calculations. 
To run an RMC simulation (with no energy constraint):
Required initial files:
-	runfile.py
-	your_epdf.gr
-	initial_structure.vasp
In this case, empirical structural constraints are necessary to prevent the Monte Carlo algorithm from producing unphysical configurations. Coordination and distance constraints are available under the Validators class. 

```
SlabThickness
Coordination: {
    BulkCoordinationRange
    SurfaceCoordinationRange
    SurfaceDistance
}
```
