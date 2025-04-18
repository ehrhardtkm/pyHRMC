# import lammps
# from lammps import LAMMPS_INT, LMP_STYLE_GLOBAL
# from lammps import PyLammps
from pymatgen.io.ase import AseAtomsAdaptor
from ase.io import write
import gc, os, psutil, subprocess
import tempfile

class Energy:
    def __init__(self, lmp_file, task_id, lmp_exec="lmp_serial"):
        """
        lmp_file: List of LAMMPS input script lines
        task_id: Unique identifier per task
        lmp_exec: Executable name or path (e.g. 'lmp_serial', 'lmp', or '/path/to/lmp')
        """
        self.lmp_file = lmp_file
        self.task_id = task_id
        self.lmp_exec = lmp_exec

    def lammps_file(self, structure, filename):
        atoms = AseAtomsAdaptor.get_atoms(structure)
        with open(filename, 'w') as f:
            write(f, atoms, format="lammps-data", specorder=structure.spec_order)

    def lammps_energy(self, structure, nsteps, success_step, dump_freq):
        with tempfile.TemporaryDirectory() as tmpdir:
            data_path = os.path.join(tmpdir, f"{self.task_id}.lmp")
            input_path = os.path.join(tmpdir, "in.lammps")

            # Write data file
            self.lammps_file(structure, data_path)

            # Write input script
            with open(input_path, "w") as f:
                for line in self.lmp_file:
                    line = line.strip()
                    if line.startswith("read_data"):
                        f.write(f"read_data {data_path}\n")
                    elif line.startswith("reset_timestep"):
                        f.write(f"variable last_step equal {success_step}\n")
                        f.write(line + "\n")
                    else:
                        f.write(line + "\n")

                # Append energy printout
                f.write("variable e equal pe\n")
                f.write("print \"${e}\"\n")  # No file output, goes to stdout

            # Run LAMMPS
            result = subprocess.run(
                [self.lmp_exec, "-in", input_path],
                cwd=tmpdir,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )

            if result.returncode != 0:
                print(f"[{self.task_id}] LAMMPS exited with code {result.returncode}")
                print(f"[{self.task_id}] STDOUT:\n{result.stdout}")
                print(f"[{self.task_id}] STDERR:\n{result.stderr}")

            # Parse stdout to extract potential energy
            pe = None
            for line in result.stdout.strip().splitlines()[-10:]:
                try:
                    pe = float(line.strip())
                    break
                except ValueError:
                    continue

            if pe is None:
                print(f"[{self.task_id}] Warning: PE not found in output. Defaulting to 0.0")
                pe = 0.0

            max_unc = 0.0  # Placeholder if you don't use MaxUnc for now

        gc.collect()
        return pe, max_unc



    # def lammps_energy(self, structure, nsteps, success_step, dump_freq):
    #     print(f"[in function before {self.task_id}] Final cleanup, open files: {len(psutil.Process(os.getpid()).open_files())}")

    #     # Step 1: Create a real temp file and write structure
    #     with tempfile.NamedTemporaryFile(mode='w', suffix='.lmp', delete=False) as tmpfile:
    #         lmp_data_path = tmpfile.name
    #         self.lammps_file(structure, lmp_data_path)

    #     pe, max_unc = 0.0, 0.0
    #     with lammps.lammps(cmdargs=["-log", "none", "-screen", "none", "-nocite"]) as lmp:
    #         if nsteps != 0:
    #             for line in self.lmp_file:
    #                 if line.startswith("read_data"):
    #                     lmp.command(f"read_data {lmp_data_path}")
    #                 elif line.startswith("reset_timestep"):
    #                     lmp.command(f"variable last_step equal {success_step}")
    #                     lmp.command(line)
    #                 elif nsteps % dump_freq == 0 and line.startswith("compute"):
    #                     max_unc = lmp.extract_compute("MaxUnc", 0, 0)
    #                 else:
    #                     lmp.command(line)

    #         elif nsteps == 0:
    #             for line in self.lmp_file:
    #                 if line.startswith("read_data"):
    #                     lmp.command(f"read_data {lmp_data_path}")
    #                 else:
    #                     lmp.command(line)

    #         lmp.command("variable e equal pe")
    #         lmp.command("run 0")

    #         energy = lmp.extract_compute("thermo_pe", LMP_STYLE_GLOBAL, LAMMPS_INT)
    #         pe = round(energy, 5)
 
    #     os.unlink(lmp_data_path)  # Unlink here as well

    #     print(f"[in function after {self.task_id}] Final cleanup, open files: {len(psutil.Process(os.getpid()).open_files())}")
        
    #     try:
    #         os.remove(lmp_data_path)
    #     except FileNotFoundError:
    #         pass

    #     gc.collect()

    #     return pe, max_unc