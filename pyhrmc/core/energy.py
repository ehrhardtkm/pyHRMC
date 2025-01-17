import lammps
from lammps import LAMMPS_INT, LMP_STYLE_GLOBAL
from lammps import PyLammps
from pymatgen.io.ase import AseAtomsAdaptor
from ase.io import write


class Energy:
    def __init__(self, lmp_file, task_id):
        self.lmp = lammps.lammps(cmdargs=["-log", "none", "-screen", "none", "-nocite"])
        self.lmp_file = lmp_file
        self.task_id = task_id

    def lammps_file(self, structure):
        atoms = AseAtomsAdaptor.get_atoms(structure)
        filename = f"{self.task_id}.lmp"
        write(filename, atoms, format="lammps-data", specorder=structure.el_list)

    def simple_lmp(self, structure):
        self.lammps_file(structure)

        for line in self.lmp_file:
            if line.startswith("read_data"):
                self.lmp.command(f"read_data {self.task_id}.lmp")
            else:
                self.lmp.command(line)
        """
        alternative is to overwrite the in.lmp copy somewhere and then read that in with lmp.file()
        """
        self.lmp.close()
        return

    def lammps_energy(self, structure, nsteps, success_step, dump_freq):
        max_unc = 0
        self.lammps_file(structure)
        if nsteps != 0:
            for line in self.lmp_file:
                if line.startswith("read_data"):
                    self.lmp.command(f"read_data {self.task_id}.lmp")
                elif line.startswith("reset_timestep"):
                    self.lmp.command(f"variable last_step equal {success_step}")
                    self.lmp.command(line)
                elif nsteps % dump_freq == 0 and line.startswith("compute"):
                    max_unc = self.lmp.extract_compute("MaxUnc", 0, 0)
                else:
                    self.lmp.command(line)

        elif nsteps == 0:
            for line in self.lmp_file:
                if line.startswith("read_data"):
                    self.lmp.command(f"read_data {self.task_id}.lmp")
                else:
                    self.lmp.command(line)

        self.lmp.command("variable e equal pe")
        self.lmp.command("run 0")

        energy = self.lmp.extract_compute("thermo_pe", LMP_STYLE_GLOBAL, LAMMPS_INT)
        pe = round(energy, 5)

        # Clean up
        self.lmp.close()
        return pe, max_unc
