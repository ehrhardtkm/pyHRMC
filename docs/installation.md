Installation
===

pyHRMC is not yet available on conda or pip. For the time being, install the package with the following steps:

1. Create and activate a virtual environment with the required packages. Perform via conda with the command:
```
conda create -n pyHRMC pip python==3.11
conda activate pyHRMC
```

2. Download a copy of the repository locally. The most recent version of the code can be found at https://github.com/ehrhardtkm/pyHRMC.
```
git clone https://github.com/ehrhardtkm/pyHRMC.git
cd pyHRMC
```

3. Install the package in development mode with pip. 
```
pip install -e .
```