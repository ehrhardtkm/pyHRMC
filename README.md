![CI](https://github.com/ehrhardtkm/pyhrmc/actions/workflows/ci.yml/badge.svg)
![License: BSD](https://img.shields.io/badge/License-BSD--3--Clause-blue.svg)
![Python Versions](https://img.shields.io/pypi/pyversions/pyhrmc)
[![Documentation](https://img.shields.io/badge/docs-pyHRMC-blue)](https://ehrhardtkm.github.io/pyHRMC/)

About 
--- 
pyHRMC is designed for HRMC simulations using experimental electron pair distribution functions as a primary constraint. This packagerelies heavily and uses code from [pymatgen](https://pymatgen.org/), which is released under the MIT license.

Full documentation can be found at https://ehrhardtkm.github.io/pyHRMC/

Installation
 --- 
Prior to installing pyHRMC, LAMMPS must be installed and built in serial. Additionally, if using a FLARE potential, LAMMPS must be compiled with FLARE. Instructions for these steps can be found at these links:

- https://docs.lammps.org/Install.html
- https://mir-group.github.io/flare/installation/lammps.html

To install pyHRMC, first create a virtual environment:
```
conda create -n pyHRMC pip python==3.11
conda activate pyHRMC
```

Installation can then be performed in the new environment. pyHRMC is currently available on PyPi for `pip install`: 
```
pip install pyhrmc
```

If users desire to modify the code from their own needs, we recommend the following steps instead:
``` 
conda create -n pyHRMC pip python==3.11 
conda activate pyHRMC
git clone https://github.com/ehrhardtkm/pyHRMC.git
cd pyHRMC
pip install -e .
```

