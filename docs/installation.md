Installation
===

Prior to installing pyHRMC, LAMMPS must be installed and built in serial. Additionally, if using a FLARE potential, LAMMPS must be compiled with FLARE. Instructions for these steps can be found at these links:

- https://docs.lammps.org/Install.html
- https://mir-group.github.io/flare/installation/lammps.html

To install pyHRMC, first create a virtual environment with the required packages. You can perform this via conda with the commands:
```
conda create -n pyHRMC pip python==3.11
conda activate pyHRMC
```

Installation can then be performed in the new environment. pyHRMC is not yet available on conda, but is available on PyPi for `pip install`: 
```
pip install pyhrmc
```

Alternatively, if users desire to modify the code from their own needs, we recommend the following steps instead:

Create an environment:
``` 
conda create -n pyHRMC pip python==3.11 
conda activate pyHRMC
```
Download a copy of the repository locally. The most recent version of the code can be found at https://github.com/ehrhardtkm/pyHRMC.
```
git clone https://github.com/ehrhardtkm/pyHRMC.git
cd pyHRMC
```


3. Install the package in development mode with pip.
```
pip install -e .
```

## Testing 
To verify that pyHRMC is correctly installed, you can run the unit tests that are available in the `tests/` directory and implemented with pytest. 

Prior to testing, ensure that you have both LAMMPS and FLARE installed to your machine. Otherwise, the HRMC test will fail. If you do not intend on using the HRMC function and are only interested in testing the RMC function, you can disregard this error and only check to see if the test succeeds when `Hybrid = False`.

From within the `pyhrmc` folder, run the following commands:
```
pip install pytest
cd tests
pytest
```