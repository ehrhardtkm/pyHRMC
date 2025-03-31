# -*- coding: utf-8 -*-

from pyhrmc.core.rmc import RMC

if __name__ == "__main__":
   
    rmc = RMC(experimental_G_csv="7ps_blur2.txt", sigma = 0.5, q_scatter = 0.9995, q_temp = 0.9995, init_temp = 2000, dump_freq = 500)
    rmc.run_rmc(
        num_processes = 16,
        initial_structure="crystal_d3_slab.vasp",
        experimental_G_csv="7ps_blur2.txt",
        keV=200,
        prdf_args={"bin_size": 0.04},
        transformations={
            "AtomHop": {}, 
        },

        validators={
            "SlabThickness": {"max_thickness": 51.5},
            "DistancesCoordination": {
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
        gaussian_blur = 2,
        charges={
            "Al": 0,
            "O": 0,
            },
        max_steps= 1000000,
    )
