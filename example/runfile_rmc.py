# -*- coding: utf-8 -*-

from pyhrmc.core.rmc_coord import RMC_coord

if __name__ == "__main__":
   
    rmc = RMC_coord(experimental_G_csv="7ps_blur2.gr", sigma = 0.01, q_scatter = 1, dump_freq = 500, hybrid = False)
    rmc.run_rmc(
        num_processes = 16,
        initial_structure="crystal_d3_slab.vasp",
        experimental_G_csv="7ps_blur2.gr",
        keV=200,
        prdf_args={"bin_size": 0.04},
        transformations={
            "AtomHop": {},
        },

        validators={
            "SlabThickness": {"max_thickness": 51.5},
            "DistancesCoordinationCoord": {
                "MinDistances": {
                    ("Al", "Al"): 2.0,
                    ("Al", "O"): 1.6,
                    ("O", "O"): 2.0
                    },
                "BulkCoordinationRange": {"Al": {"Al" : [0, 0], "O": [4, 5]}, "O": {"Al": [3, 4], "O": [0, 0]} },
                "SurfaceCoordinationRange": {"Al": {"Al": [0, 0], "O": [4, 5]}, "O": {"Al": [2, 4], "O": [0, 0]} },
                "SurfaceDistance": 3
                    }
            },
        qmin = 2,
        charges={
            "Al": 0,
            "O": 0,
            },
        max_steps= 1000000,
    )
