pyHRMC.validators
===
Validators are used to impose physical constraints on the Monte Carlo algorithm and enforce the evolution of physically realistic structural configurations. In traditional RMC, these constraints are the sole avenue in imposing "chemical knowledge" on the simulation. In HRMC, the need for accurate and restrictive constraints are alleviated because the structure's energy is also used to determine whether to accept or reject a proposed step.

Since no set of validators is universally applicable across systems, users must select values that are appropriate for their use case. In the event of performing a bulk HRMC simulation instead of a slab, ...*EDIT CODE TO PERMIT BULK SIMULATIONS*

```
SlabThickness
Coordination: {
    BulkCoordinationRange
    SurfaceCoordinationRange
    SurfaceDistance
}
```
