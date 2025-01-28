Bulk vs slab simulations
===
pyHRMC is designed to handle bulk cells or slab cells in a vacuum. While the package will automatically detect the cell type and handle periodic boundary conditions, several considerations and limitations exist. 

For a slab cell, pyHRMC assumes that the vacuum space above and below the cell exist in the z-direction, and that the structure extends fully to the cell boundaries in the x and y directions. Additionally, the slab must be orthogonal to the z-axis for the package to correctly detect the cell thickness. When running a slab cell simulation and employing the `Validators.DistancesCoordination` constraint, be sure to set constraints for `SurfaceCoordinationRange` and `SurfaceDistance` in addition to `BulkCoordinationRange`. 

For example:
```
"DistancesCoordination": {
    "MinDistances": {
        ("Al", "Al"): 2.0,
        ("Al", "O"): 1.6,
        ("O", "O"): 2.0
        },
    "BulkCoordinationRange": {"Al": {"Al" : [0, 0], "O": [4, 6]}, "O": {"Al": [3, 4], "O": [0, 0]} },
    "SurfaceCoordinationRange": {"Al": {"Al": [0, 0], "O": [4, 6]}, "O": {"Al": [2, 4], "O": [0, 0]} },
    "SurfaceDistance": 3
        }

```

In the event of performing a bulk HRMC simulation instead of a slab, the cell must also be orthogonal in the z-axis and assumes that there is **not** vacuum space in the z-direction. If the difference between the lattice dimension in the z-direction and the structure thickness along z is greater than 1, the package will incorrectly assume that the cell is a slab and may perform erroneously. The `Validators.DistancesCoordination` constraint should be used with only `BulkCoordinationRange`.