from setuptools import find_packages, setup

setup(
    name="pyHRMC",
    version="0.0.0",
    description="Tools to run HRMC with experimental ePDF data using Simmate",
    url="https://github.com/ehrhardtkm/pyHRMC",
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
	'simmate'
    ],
    author="Karen M. Ehrhardt",
    author_email="kme@unc.email.edu"
)
