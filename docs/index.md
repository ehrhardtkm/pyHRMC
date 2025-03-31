# Welcome to pyHRMC!

pyHRMC is a tool designed to aid users in perform Hybrid Reverse Monte Carlo simulations from electron-based scattering. Some key features of the package include:

1. Flexible modeling capabilities, including simulation of slab structures
2. Accurate scattering cross sections for electron scattering
3. Machine-learned potentials to predict structural energies with high accuracy 
4. Highly parallelizable architecture for major speedups  

## Getting started
This site provides the official documentation for pyHRMC. 

To install pyHRMC and start your first simulation, please see the [Installation](https://ehrhardtkm.github.io/pyHRMC/installation/) steps.

## About HRMC

Hybrid Reverse Monte Carlo (HRMC) is a computational method that stochastically perturbs the atomic configuration of a structure to fit its configuration to experimental structural data. Typically, this data is a pair distribution function (PDF) that is collected from x-ray, neutron, or electron sources. In particular, electron scattering is advantageous for nanoscale and ultrathin materials because of the high spatial resolution and larger scattering cross-section of electrons. As materials become increasingly nanoscale, the ability to experimentally collect scattering information from disordered materials is useful. However, PDFs provide limited information on the atomistic arrangement of atoms.

This is where Reverse Monte Carlo (RMC) comes in. In traditional implementations of RMC, the simulation is designed to gradually fit a computed PDF of a cell to an experimental PDF by way of random atomic hops, according to the standard Metropolis Monte Carlo algorithm. While this approach is effective at matching the experimental and computed PDFs, it can fail to produce physically realistic structures. This flaw is due to the stochastic nature of the algorithm, and the limited chemical information that guides the simulation. Often, users specify coordination number and bond angle ranges, minimum bond distances, and more to try to constrain the simulation sufficiently to produce realistic structures. However, these inputs require the user to possess a high degree of understanding about the system that they are trying to model and the accuracy of the simulated structure depends heavily on these parameters. 

To mitigate the issue of user-selected constraint parameters, HRMC was introduced. In addition to fitting the compute PDF to experimental data, HRMC also includes a Boltzmann-type probability term to account for the energy of the structure. This approach therefore effectively minimizes the structural energy while finding and experimental fit, resulting in more flexible constraints to better allow transition states during atom hopping and has been shown to more accurately reproduce chemically accurate structures.

pyHRMC is a Python implementation of this HRMC algorithm, with several added features. 


## Resources

This code was developed with tools from [pymatgen](https://pymatgen.org/) and [matminer](https://hackingmaterials.lbl.gov/matminer/). 

For bugs and other issues, please submit them on [Github](https://github.com/ehrhardtkm/pyHRMC).

Questions or concerns? Contact kme@email.unc.edu

### Citing pyHRMC
We are currently preparing a publication. This section will be updated when the publication is released.

## Compatibility
pyHRMC is currently comptability with Python <= 3.11. 

## License
pyHRMC is licensed under a 3-clause BSD license. Portions of the code were borrowed from pymatgen, matminer, and Simmate, which are also licensed under MIT and BSD-style licenses. For details, please see [LICENSE](https://github.com/ehrhardtkm/pyHRMC/LICENSE.md).


