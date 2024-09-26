Non-hybrid simulations (RMC)
===
Traditional RMC simulations can also be performed with pyHRMC. It’s important to note that traditional RMC is typically less accurate than HRMC, and for this reason, we don’t recommend it. However, it can be useful in some cases, such as when a suitable force field is not available for energy calculations. In this case, empirical structural constraints are necessary to prevent the Monte Carlo algorithm from producing unphysical configurations. Coordination and distance constraints are available under the Validators class. 

To run an RMC simulation (with no energy constraint), use the following input files:

- runfile.py
- your_epdf.gr
- initial_structure.vasp

The runfile must be configured to appropriately to avoid running features designing for HRMC. When instantiating the RMC class, set `hybrid=False` (defaults to True). Additionally, `q_temp` and `init_temp` arguments are not required for non-hybrid RMC (defaults to None):

Below is an example of a non-hybrid RMC setup with no quenching of the error constraint parameter `sigma`:
```
    rmc = RMC(
        experimental_G_csv="al2o3_5nm_gr.txt", 
        sigma = 0.05, 
        q_scatter = 1, 
        hybrid= False
        )
```

When calling the `run_rmc` method, the configuration is the same as HRMC. However, we strongly suggest using additional validators to appropriately constrain the system.