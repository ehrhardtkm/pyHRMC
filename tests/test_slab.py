
import os
import pytest
from pyhrmc.core.rmc import RMC


# Define the path to the data folder
DATA_FOLDER = os.path.join(os.path.dirname(__file__), 'data')

# Define the path to the data folder
TEST_FOLDER = os.path.join(os.path.dirname(__file__), 'tests')

# Directly specify the two files you want to use
input_files = {
    'al2o3_5nm_gr.txt': os.path.join(DATA_FOLDER, 'al2o3_5nm_gr.txt'),
    'POSCAR_47A-d33': os.path.join(DATA_FOLDER, 'POSCAR_47A-d33'),
    'in.lmp': os.path.join(DATA_FOLDER, 'in.lmp'),
    'in_init.lmp': os.path.join(DATA_FOLDER, 'in_init.lmp'),
    'in_accept.lmp': os.path.join(DATA_FOLDER, 'in_accept.lmp')
}


# Cleanup fixture
@pytest.fixture
def cleanup_output_files():
    yield
    # Code after yield runs after the test  
    test_file = "test_slab.py"

    print("Current working directory:", os.getcwd())
    # After the test, remove all files in the current directory except the test file
    for file in TEST_FOLDER:
        if os.path.isfile(file) and file != test_file:
            try:
                os.remove(file)
            except Exception as e:
                print(f"Error removing file {file}: {e}")


@pytest.mark.parametrize("hybrid", [True, False], )
def test_hrmc_std(hybrid, capfd, cleanup_output_files):
    rmc = RMC(
        experimental_G_csv=input_files['al2o3_5nm_gr.txt'], 
        sigma = 0.05, 
        q_scatter = 1, 
        q_temp = 0.99995, 
        init_temp = 1000, 
        hybrid = hybrid
        )
    
    max_steps = 10
    rmc.run_rmc(
        num_processes = 1,
        initial_structure=input_files['POSCAR_47A-d33'],
        keV=200,
        prdf_args={"bin_size": 0.04},
        # error_constant=0.05,
        transformations={
            # "ASECoordinatePerturbation":{}, #simmate version of RattleMutation from ase, rattle_prob auto set to 0.3
            "AtomHop": {}, # consider adding a second, smaller step size "max_step": 0.2
        },

        validators={
            "SlabThickness": {"max_thickness": 47.4},
            "DistancesCoordination": {
                "MinDistances": {
                    ("Al", "Al"): 1.6,
                    ("Al", "O"): 1.6,
                    ("O", "O"): 1.6
                },
                "BulkCoordinationRange": {"Al": {"Al" : [0, 2], "O": [2, 7]}, "O": {"Al": [1, 4], "O": [0, 1]} },
                "SurfaceCoordinationRange": {"Al": {"Al": [0, 2], "O": [1, 7]}, "O": {"Al": [1, 4], "O": [0, 1]} },
                "SurfaceDistance": 5
                }
            },

        lmp_init=input_files['in_init.lmp'],
        lmp_test = input_files['in.lmp'],
        lmp_accept=input_files['in_accept.lmp'],
        max_steps= max_steps
    )

    # Now, read the contents of one of the output files to verify its correctness
    output_file = "error_plotting.txt"
    
    # Ensure the file exists and read its content
    assert os.path.exists(output_file), f"{output_file} file does not exist"

    # Capture printed output using capfd
    captured = capfd.readouterr()
    # Validate the content of the file (example)
    expected1 = "counter"
    assert expected1 in captured.out
    expected2 = "moves"
    assert expected2 in captured.out

    pass