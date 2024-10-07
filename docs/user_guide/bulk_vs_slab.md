Bulk vs slab simulations
===
pyHRMC is designed to handle bulk cells or slab cells in a vacuum. While the package will automatically detect the cell type and handle periodic boundary conditions, several considerations and limitations exist. 

For a slab cell, pyHRMC assumes that the vacuum space above and below the cell exist in the z-direction, and that the structure extends fully to the cell boundaries in the x and y directions. Additionally, the slab must be orthogonal to the z-axis for the package to correctly detect the cell thickness. When running a slab cell simulation and employing the `Validators.coordination` constraint, be sure to set constraints for `SurfaceCoordinationRange` and `SurfaceDistance` in addition to `BulkCoordinationRange`. 

For example:
```
"Coordination": {
    "BulkCoordinationRange": {
        "Al": {"Al" : [0, 2], "O": [2, 7]}, 
        "O": {"Al": [1, 4], "O": [0, 1]} 
        },
    "SurfaceCoordinationRange": {
        "Al": {"Al": [0, 2], "O": [1, 7]}, 
        "O": {"Al": [1, 4], "O": [0, 1]} 
        },
    "SurfaceDistance": 5
    }
```

In the event of performing a bulk HRMC simulation instead of a slab, the cell must also be orthogonal in the z-axis and assumes that there is **not** vacuum space in the z-direction. If the difference between the lattice dimension in the z-direction and the structure thickness along z is greater than 1, the package will incorrectly assume that the cell is a slab and may perform erroneously. The `Validators.coordination` constraint should be used with only `BulkCoordinationRange`.