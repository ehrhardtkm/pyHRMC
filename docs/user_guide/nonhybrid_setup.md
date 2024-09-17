Non-hybrid simulations (RMC)
===
Traditional RMC simulations can also be performed with pyHRMC. It’s important to note that traditional RMC is typically less accurate than HRMC, and for this reason, we don’t recommend it. However, it can be useful in some cases, such as when a suitable force field is not available for energy calculations. 
To run an RMC simulation (with no energy constraint), use the following input files:

- runfile.py
- your_epdf.gr
- initial_structure.vasp

In this case, empirical structural constraints are necessary to prevent the Monte Carlo algorithm from producing unphysical configurations. Coordination and distance constraints are available under the Validators class. 

