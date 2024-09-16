from pymatgen.io.lammps.data import LammpsData
import multiprocessing
import lammps
from lammps import LAMMPS_INT, LMP_STYLE_GLOBAL


class Lammps_HRMC():
    @staticmethod
    def lammps_file(structure):
        lammps = LammpsData.from_structure(structure)
        lammps.set_charge_atom_type({'Al':1.4175, 'O':-0.945})
        lammps.write_file('al2o3_lammps.lmp')
         
    
    def lammps_energy(self, structure):
        #Create lammps file
        self.lammps_file(structure)
        
        # Initialize LAMMPS instance
        lmp = lammps.lammps(cmdargs=["-log", "none", "-screen", "none",  "-nocite"])
        
        # # lines = open('in.lmp','r').readlines()
        # # for line in lines: lmp.command(line)
        lmp.file('in.lmp')
        lmp.command("variable e equal pe")
        
        # run 0 to get energy of perfect lattice
        # emin = minimum energy
        
        lmp.command("run 0")
        
        natoms = lmp.extract_global("natoms")
        energy = lmp.extract_compute("thermo_pe",LMP_STYLE_GLOBAL,LAMMPS_INT) / natoms
        
        pe = round(energy, 5) 

        # Clean up
        lmp.close() 
        return pe


