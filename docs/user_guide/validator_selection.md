pyHRMC.validators
===
Validators are used to impose physical constraints on the Monte Carlo algorithm and enforce the evolution of physically realistic structural configurations. In traditional RMC, these constraints are the sole avenue in imposing "chemical knowledge" on the simulation. In HRMC, the need for accurate and restrictive constraints are alleviated because the structure's energy is also used to determine whether to accept or reject a proposed step.

Since no set of validators is universally applicable across systems, users must select values that are appropriate for their use case. 

Validators.SlabThickness
---
```
"SlabThickness": {"max_thickness": float}
```
When setting a `SlabThickness` constraint, a `max_thickness` value must be set. This is the max thickness that the simulation will allow a slab cell to reach in the z-direction (note that this is thickness, not the max Cartesian coordinate). This constraint should not be used when running bulk cell simulations, since cell length along z should be roughly equal to the thickness.

Validators.Coordination
---
```
"Coordination": {
    "BulkCoordinationRange",
    "SurfaceCoordinationRange" = None,
    "SurfaceDistance" = None
}
```

For the `BulkCoordinationRange` and `SurfaceCoordinationRange`, the `Coordination` constraint expects a dictionary for each species type. Within each dictionary, a subdictionary with the minimum and maximum allowable coordination numbers should be included for each species type in the structure. This should look like a nested dictionary. Below is an example for alumina (Al and O species):
```
"BulkCoordinationRange": {
    "Al": {"Al" : [0, 2], "O": [2, 7]}, 
    "O": {"Al": [1, 4], "O": [0, 1]} 
    }
``` 
The `SurfaceDistance` can be an integer or float and is used to differentiate surface and bulk atoms. This value is the distance from the top and bottom of the slab cell which contains the region of atoms that are constrained by `SurfaceCoordiantionRange`. All remaining atoms between these two surface regions are constrained by `BulkCoordinationRange`.

Validators.SiteDistance
---
```
"SiteDistance": {"distance_cutoff": }
```

Validators.TargetDensity
---
```
"TargetDensity": {        
    target_density,
    percent_allowance=0.05
}
```


Example
---
An example of a `run_rmc` call with the `SlabThickness` and 'Coordination` validators in included in the documentation for [runfile setup](https://ehrhardtkm.github.io/pyHRMC/user_guide/getting_started/runfile_setup/).