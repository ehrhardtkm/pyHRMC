Outputs
===

Once the HRMC simulation is started, pyHRMC will generate several files that it will continuously update during the simulation. These files contain data regarding the simulation and can be used for post-processing:

## error_plotting.txt

With each acceptance of a batch of steps, the updated information will be appended to the `error_plotting.txt` file:

```
{Iteration_number} {Error} {Energy_per_atom} {Sigma} {Temperature} {Maximum_uncertainty}
```

- **Iteration number:**  the number of HRMC loops performed. To obtain the total number of steps, multiple this value by the number of processes used in each iteration
- **Error:** the MSE of the experimental and calculated ePDF
- **Energy per atom:** the units of energy will depend on the units specified in the LAMMPS input script.
- **Sigma:** this value is used to calculate the probability of accepting a step that increases the error (i.e. is the worse fit to experimental ePDF). Higher values correspond to an increased likelihood of acceptance, while a value of 0 means that such a step is never accepted. If a value of q_scatter < 1 is used, this value will gradually decrease over the course of the simulation. For details, see [Background and Theory](https://ehrhardtkm.github.io/pyHRMC/background_theory/)
- **Temperature:** this value is used to calculate the Boltzmann-type probability of accepting an energetically unfavorable step. Higher temperatures correspond to an increased probability of acceptance. If a value of q_temp < 1 is used, this value will gradually decrease over the course of the simulation. For details, see [Background and Theory](https://ehrhardtkm.github.io/pyHRMC/background_theory/)
- **Max uncertainty:** this is an **optional** output which can be set by including the appropriate lines in the [LAMMPS input script](https://ehrhardtkm.github.io/pyHRMC/user_guide/getting_started/lammps_input). If turned on, the maximum uncertainty computed from the FLARE potential will be printed according to the `dump_freq` that is set during the creation of the RMC object. Otherwise, this column will be populated by zeroes. 

In the event of a non-hybrid simulation, `error_plotting.txt` will include the following information instead:
```
{Iteration_number} {Error} {Sigma} 
```

## output.vasp

The most current structure in the POSCAR format. This file is overwritten every time an [iteration](https://ehrhardtkm.github.io/pyHRMC/runfile_setup/) of the HRMC loop is accepted.

## XDATCAR:

A trajectory file of the simulation in the XDATCAR format. Dumped every iteration of the HRMC loop, and therefore every number of steps equal to the number of processes used.

## pdf.txt

A text file of the calculated ePDF, which is overwritten with a frequency of `dump_freq` iterations. The file includes the radius (x) values followed by the G(r) (y) values.

## pdfs.png

A plotted figure of the experimental vs. calculated ePDFs. This is plotted using the same data that is included in the `pdf.txt` file and overwritten at the same frequency, but is saved as a figure for easy viewing purposes.

***
## Other files

Additionally, several files are created and used by HRMC during the simulation, but do not provide useful information to the user:

- {int}.lmp

- last_step.txt

This files are used to maintain continuity between calls to LAMMPS and to distinguish the structures between processes.
