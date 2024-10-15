pyHRMC.validators
===
Validators are used to impose physical constraints on the Monte Carlo algorithm and enforce the evolution of physically realistic structural configurations. In traditional RMC, these constraints are the sole avenue in imposing "chemical knowledge" on the simulation. In HRMC, the need for accurate and restrictive constraints are alleviated because the structure's energy is also used to determine whether to accept or reject a proposed step.

Since no set of validators is universally applicable across systems, users must select values that are appropriate for their use case. 

SlabThickness
---
```
"SlabThickness": {"max_thickness": float}
```
When setting a `SlabThickness` constraint, a `max_thickness` value must be set. This is the max thickness that the simulation will allow a slab cell to reach in the z-direction (note that this is thickness, not the max Cartesian coordinate). This constraint should not be used when running bulk cell simulations, since cell length along z should be roughly equal to the thickness.

Coordination
---
```
"Coordination": {
    "BulkCoordinationRange": {
        "species": {
            "species1" : [int(min_coordination), int(max_coordination)], 
            "species2": [int(min_coordination), int(max_coordination)], 
            ...
            }
        ...
        },
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

DistancesCoordination
---
```
"DistancesCoordination": {
    "MinDistances": {
        ("species", "species"): float,
        ...
    },
    "BulkCoordinationRange": {
        "species": {
            "species1" : [int(min_coordination), int(max_coordination)], 
            "species2": [int(min_coordination), int(max_coordination)], 
            ...
            }
        ...
        },
    "SurfaceCoordinationRange" = None
    "SurfaceDistance"= None
    }
```

In addition to the same coordination check that is executed by the `Coordination` constraint, `DistanceCoordination` will conduct a check to determine if any interatomic distances are below the user-defined minimum distances. Distance cutoffs for this constrain should be defined pairwise for the species types. If any distances below the corresponding cutoff is found, the step will automatically be rejected.

When checking the MinDistances, the program will recognize the pairwise combination regardless of order (e.g. ("Al", "O") is the same as ("O", "Al")). Below is an example for alumina (Al and O species):
```
"DistancesCoordination": {
    "MinDistances": {
        ("Al", "Al"): 1.6,
        ("Al", "O"): 1.6,
        ("O", "O"): 1.8
    },
    "BulkCoordinationRange": {
        "Al": {"Al" : [0, 2], "O": [2, 7]}, 
        "O": {"Al": [1, 4], "O": [0, 1]} 
        },
    "SurfaceCoordinationRange": {
        "Al": {"Al": [0, 2], "O": [1, 7]},
        "O": {"Al": [1, 4], "O": [0, 1]} 
        },
    "SurfaceDistance": 3
    }
```

If the same minimum distance can be applied to all element pairs, the `SiteDistance` constraint will be faster and can be used in conjunction with the `Coordination` constraint. 


SiteDistance
---
```
"SiteDistance": {
    distance_cutoff = float
}
```
The `SiteDistance` constraint will uses the [distance_matrix](https://pymatgen.org/pymatgen.core.html) property from pymatgen to check all interatomic distances in the system against a single, universal minimum distance cutoff. If any interatomic distances are found to be less than the user-specified cutoff, the program will automatically reject the step.

TargetDensity
---
```
"TargetDensity": {        
    "target_density": float,
    "percent_allowance": 0.05,
    "check_slab": True
}
```
The structure density is compared to the user-defined target density. If the density does not fall within the percent allowance of the target density (defaults to 5%), the step will automatically be rejected. If using a bulk cell, use `"check_slab": False`, otherwise this will default to True.


Example
---
An example of a `run_rmc` call with the `SlabThickness` and 'Coordination` validators in included in the documentation for [runfile setup](https://ehrhardtkm.github.io/pyHRMC/user_guide/getting_started/runfile_setup/).