# -*- coding: utf-8 -*-
import os
import numpy as np
from numpy import exp
from random import random
import pyHRMC.transformers as transform_mod
import pyHRMC.validators as validator_mod
from pyHRMC.core.slab import RdfSlab
from pyHRMC.core.interpolator import CrossSection
from pyHRMC.core.hrmc import Lammps_HRMC

import multiprocessing
from multiprocessing.shared_memory import SharedMemory
import copy
from itertools import combinations
import pickle
import warnings
from pymatgen.core import Structure

class RMC():
      
    def __init__(self, experimental_G_csv, sigma, q_scatter, q_temp = None, init_temp=None, hybrid=True):
        self.hybrid = hybrid
        self.experimental_G_csv = experimental_G_csv
        self.batched_error_constant = sigma
        self.batched_temp = init_temp
        self.current_error = float()
        self.current_energy = float()
        self.q_temp = q_temp
        self.q_scatter = q_scatter
        self.nsteps = 0
        self.success_step = 0

        if self.hybrid == True and self.q_temp == None or self.batched_temp == None:
            raise RuntimeError(
                "If running HRMC, please provide values for q_temp and init_temp.\n"
                "If running RMC, specify that hybrid == False\n"
            )
    
    def apply_oxi_state(self, valences, struct):
        # BV = BVAnalyzer()
        # valences = BV.get_valences(self)
        for i in range(len(struct.sites)):
            struct.sites[i].oxi_state = valences[i]
        return
    
    
    #generate and check function to parallelize
    def make_candidate_structure(self, current_structure, validator_objects, transformer, task_id, lmp_input=None):
        new_structure = transformer.apply_transformation(current_structure)
        new_structure.xyz_df = new_structure.xyz()
    
        is_valid = True  # true until proven otherwise
        for validator in validator_objects:
            if not validator.check_structure(new_structure):
                is_valid = False
                break  # one failure is enough to stop and reject the structure
        if not is_valid:
            new_structure = False
        else:
            #  do pdf error check
            neighborlist = new_structure.get_all_neighbors(r=10.0)
            if self.hybrid == True:
                keep_new, new_error, max_unc = self.choose_acceptance("single", new_structure, neighborlist,self.batched_error_constant, self.batched_temp, lmp_input, task_id)
            else:
                keep_new = self.choose_acceptance("single", new_structure, neighborlist,self.batched_error_constant)
            if keep_new:
                pass
            else:
                new_structure = False
        return new_structure
      
    def choose_acceptance(self, version, structure, neighborlist, error_constant, temp=None, lmp_input=None, task_id=None):
        max_unc = 0
        if self.hybrid == True:
            #HRMC
            lammps_run = Lammps_HRMC(lmp_input, task_id)
            new_energy, max_unc = lammps_run.lammps_energy(structure, self.nsteps, self.success_step)
        
        structure.load_experimental_from_file(self.experimental_G_csv)
        new_error,slope = structure.prdf_error(neighborlist)
        
        keep_new = True

        """
        REORGANIZE TO BE BASED ON IF HBRID == TRUE
        """



        if new_error < self.current_error: 
            if self.hybrid == True:
                if new_energy <= self.current_energy:
                    # accept new structure
                    with open('acceptance_probabilities.txt', 'a') as out:
                        out.write(f'{version}: {keep_new}, n/a, {new_error}, {new_energy}\n')
                    pass
                elif new_energy > self.current_energy:
                    # HRMC, with energy term
                    new_h_error = ( new_error / (error_constant ** 2) ) +  ( new_energy / ((1.987204259 * 10**-3) * temp) )
                    old_h_error = ( self.current_error / (error_constant ** 2) ) +  ( self.current_energy / ((1.987204259 * 10**-3) * temp) )
                    h_prob = exp(old_h_error - new_h_error)
                    if random() < h_prob:
                        pass
                    else:
                        keep_new = False
                    #logging probability info for record keeping
                    with open('acceptance_probabilities.txt', 'a') as out:
                        out.write(f'{version}: {keep_new}, {h_prob}, {new_error}, {new_energy}\n')
                return keep_new, new_error, max_unc
            #RMC
            elif self.hybrid == False:
                # accept new structure
                with open('acceptance_probabilities.txt', 'a') as out:
                    out.write(f'{version}: {keep_new}, n/a, {new_error}\n')
                return keep_new, new_error

        else:
            if self.hybrid == True:
                # HRMC, with energy term
                new_h_error = ( new_error / (error_constant ** 2) ) +  ( new_energy / ((1.987204259 * 10**-3) * temp) )
                old_h_error = ( self.current_error / (error_constant ** 2) ) +  ( self.current_energy / ((1.987204259 * 10**-3) * temp) )
                h_prob = exp(old_h_error - new_h_error)
                if random() < h_prob:
                    pass
                else:
                    keep_new = False
                #logging probability info for record keeping
                with open('acceptance_probabilities.txt', 'a') as out:
                    out.write(f'{version}: {keep_new}, {h_prob}, {new_error}, {new_energy}\n')
                return keep_new, new_error, max_unc
            elif self.hybrid == False:
                probability = exp(-1 * (abs(new_error - self.current_error)) / error_constant)
                if random() < probability:
                    pass
                else:
                    keep_new = False
                return keep_new, new_error
      
        # neighborlist update with custom approach. then do pdf error check. use __init__ for experimental_G_csv. Then use probability and return true or False 
    def update_neighbors(self, shm_name, size, n_idx):
        #put in the parallelized method, deserialization
        #attach to existing shared memory
        shm = SharedMemory(name=shm_name)
        structure_data = shm.buf[:size].tobytes()
        structure = pickle.loads(structure_data)
    
        n_site = structure.sites[n_idx]
        nns = structure.get_sites_in_sphere(n_site.coords, r = 10.0)
        #remove the site itself because the pymatgen method includes it
        updated_neighbors = (n_idx, [nn for nn in nns if nn.index != n_idx])
        
        shm.close()
    
        return updated_neighbors      
          
    def worker_task(self, structure, lmp_input, task_id):
        lammps_run = Lammps_HRMC(lmp_input, task_id)
        energy, max_unc = lammps_run.lammps_energy(structure, self.nsteps, self.success_step)
        return energy, max_unc     
          

    def store_symbol(self, structure):
        el_switch = []
        j=0
        el_list= [structure.sites[0].specie.symbol]
        for i in range(len(structure.sites)):
            if  structure.atomic_numbers[i] != structure.atomic_numbers[j]:
                el_switch.append(i)
                el_list.append(structure.sites[i].specie.symbol)
            j=i
        return el_switch, el_list
    
    def run_rmc(self,
        num_processes = int(),
        initial_structure="amorphous_al2o3.vasp",
        experimental_G_csv="al2o3_5nm_gr.txt",
        keV=200,
        prdf_args={"bin_size": 0.04},
        #for testing, changed error_constant from 5E-7 to 0.002
        # self.batched_error_constant=0.002,
        transformations: dict = {
            # "ASECoordinatePerturbation":{}, #simmate version of RattleMutation from ase, rattle_prob auto set to 0.3
            "AtomHop": {}, # consider adding a second, smaller step size "max_step": 0.2
        },
        validators={
        },
        lmp_init='in_init.lmp',
        lmp_test = 'in.lmp',
        lmp_accept='in_accept.lmp',
        
        
        max_steps= 100,
    ):

        if os.path.exists("error_plotting.txt"):
            if self.hybrid == True:
                try:
                    with open("last_step.txt", 'r') as file:
                        restart_success_step = int(file.read())
                        self.success_step = restart_success_step
                except IndexError:
                    raise RuntimeError(
                        "You're likely trying to restart a simulation and switch between hybrid and non-hybrid RMC."
                        "Please change your runfile.py settings for the same type of simulation as previously."
                    )

                with open('error_plotting.txt', 'r') as file:
                    lines = file.readlines()
                    last_line = lines[-1].strip()
                    variables = last_line.split()
                    self.nsteps = int(variables[0])
                    try:
                        self.batched_error_constant = float(variables[3])
                        self.batched_temp = float(variables[4])
                    except IndexError:
                        raise RuntimeError(
                            "You're likely trying to restart a simulation and switch between hybrid and non-hybrid RMC.\n"
                            "Please change your runfile.py settings for the same type of simulation as previously."
                        )
            elif self.hybrid == False:
                with open('error_plotting.txt', 'r') as file:
                    lines = file.readlines()
                    last_line = lines[-1].strip()
                    variables = last_line.split()
                    self.nsteps = int(variables[0])
                    self.batched_error_constant = float(variables[2])

            warnings.warn(
                "\nA previous run has been detected. The simulation will resume at:\n"
                f"{self.nsteps} steps\n"
                "pyHRMC will automatically append the output of this run to existing output files.\n"
                "If this is undesired, please remove the existing file and restart the simulation.")

            if not os.path.exists("errors.txt"):
                warnings.warn(
                    "errors.txt is missing.")
            if os.path.exists("XDATCAR"):
                warnings.warn(
                    "XDATCAR is missing.")
            if os.path.exists("acceptance_probabilities.txt"):
                warnings.warn(
                    "acceptance_probabilities is missing.")
            


        
        """
        STAGING
        """
        
        #energy weighting staging
        final_temp = 100

        
        #batch number staging
        num_batch = num_processes
        # decay_const = 0.1
        # batch_decay = exp(- decay_const * num_batch)

        """
        INITIAL STRUCTURE AND CALCULATION SETUP
        """

        if self.hybrid == True:
            #lammps input setup
            lmp_init_in = []    
            with open(lmp_init, 'r') as file:
                for line in file:
                    lmp_init_in.append(line)
                    
            lmp_input = []    
            with open(lmp_test, 'r') as file:
                for line in file:
                    lmp_input.append(line)
                    
            lmp_acc_in = []    
            with open(lmp_accept, 'r') as file:
                for line in file:
                    lmp_acc_in.append(line)
                
                
        # convert to our python class
        initial_structure = RdfSlab.from_file(initial_structure)
        initial_structure.xyz_df = initial_structure.xyz()
        #create neighborlist
        initial_structure_neighborlist = initial_structure.get_all_neighbors(r=10.0)
        num_atoms = len(initial_structure.sites)
        
        #cache/store interpolated radii
        struct_consts = CrossSection(initial_structure)

        valences = initial_structure.oxidation_state_list()
        self.apply_oxi_state(valences, initial_structure)
        charges = struct_consts.partial_charges()
        database_radii = struct_consts.database_ionic_radii()
        interpolated_radii = struct_consts.interpolated_ionic_radii(charges, database_radii)
        setattr(initial_structure, 'interpolated_radii', interpolated_radii)

        TCS = {}
        el_tuple = initial_structure.symbol_set
        for el in el_tuple:
            el_charge = charges[el]
            TCS[el] = struct_consts.interpolated_TCS(el, el_charge, keV)
        
        #BUG: FOR TROUBLEHSOOTING ONLY, FIX THIS 
        # setattr(initial_structure, 'TCSs', TCS)
        CROSS_SECTIONS_CONSTANTS = {'Al': 5.704193e-2, 'O': 2.029453e-2}
        setattr(initial_structure, 'TCSs', CROSS_SECTIONS_CONSTANTS)
        print(initial_structure.TCSs)

    
        # load experimental G(r)
        initial_structure.load_experimental_from_file(self.experimental_G_csv)
    
        # convert validators + transformations all to proper python class
        # below is an example loop with transformations
        transformation_objects = []
        for t_name, t_kwargs in transformations.items():
            t_class = getattr(transform_mod, t_name)
            t_obj = t_class(**t_kwargs)
            transformation_objects.append(t_obj)
            
    
        # dynamically build our validator objects
        validator_objects = []
        for v_name, v_kwargs in validators.items():
            if hasattr(validator_mod, v_name):
                v_class = getattr(validator_mod, v_name)
                v_obj = v_class(**v_kwargs)
            else:
                raise Exception(f"Unknown validator provided ({v_name})")
            validator_objects.append(v_obj)

        
        """
        INITIALIZE RMC LOOP
        """
        current_structure = copy.deepcopy(initial_structure)
        el_switch, el_list = self.store_symbol(current_structure)
        current_structure.el_switch = el_switch
        current_structure.el_list = el_list
        current_structure_neighborlist = initial_structure_neighborlist
        self.current_error,slope = current_structure.prdf_error(current_structure_neighborlist)
        current_e = round(self.current_error, 5)
        #additional steps to initialize HRMC
        if self.hybrid == True:
            default_id = 100
            initial_e, max_unc = self.worker_task(initial_structure, lmp_init_in, default_id)
            self.current_energy = initial_e
            print(f'Step {self.nsteps}. Sum of residuals = {current_e}. Energy = {self.current_energy}.')
        else:
            print(f'Step {self.nsteps}. Sum of residuals = {current_e}.')
        
        counter = []
        error_list = []
        moves = 0
        moves_attempted = 0

        current_structure.plot_pdf(current_structure_neighborlist, experimental_G_csv, slope)
        
        """
        RUN RMC LOOP
        """
        while self.nsteps < max_steps:
            self.nsteps += 1
            if self.nsteps % 5000 == 0:
                if os.path.exists("pdfs.png"):
                    os.remove("pdfs.png") 
                current_structure.plot_pdf(current_structure_neighborlist, experimental_G_csv, slope)


            transformer = transformation_objects[0]  # just grab the first transformation
            
            trial_structure = copy.deepcopy(current_structure)
            new_structure = copy.deepcopy(current_structure)
            
            with multiprocessing.Pool(processes=num_processes) as pool:
                tasks = list(range(num_batch))
                if self.hybrid == True:
                    inputs = [(trial_structure, validator_objects, transformer, task_id, lmp_input) for task_id in tasks]
                else:
                    inputs = [(trial_structure, validator_objects, transformer, task_id) for task_id in tasks]
                results = pool.starmap(self.make_candidate_structure, inputs)
                                
            moves_attempted += num_batch       
            moved_atoms = []

            for candidate_structure in results:
                if candidate_structure:
                    # new_energy = energy
                    idx = candidate_structure.move_indices[0]
                    moved_atoms.append([idx, candidate_structure.sites[idx]])
                    
            # check distance between this atom and all other moved atoms
            for i, j in combinations(moved_atoms, 2):
                # distance = np.linalg.norm(i[2]-j[2])
                #this becomes:
                distance = np.linalg.norm(i[1].coords-j[1].coords)
                try: #only do this if min distance constraint is in use
                    min_distance = validator_objects[1].min_distances[f'{i[1].species_string}', f'{j[1].species_string}']
                    if distance < min_distance:
                        try:
                            moved_atoms.remove(i)
                        except:
                            pass
                except: 
                    pass
            if len(moved_atoms) == 0:
                continue

            all_changed = set()
            for c_idx, c_site in moved_atoms:
                new_structure.sites[c_idx] = c_site
            #edit neighbor list here
            #BUG: switch to individuvudal neighbor method for speedup
            new_structure_neighborlist = new_structure.get_all_neighbors(r=10.0)
            #update xyz_df
            new_structure.xyz_df = new_structure.xyz()
            
            if self.hybrid == True:
                keep_new, new_error, max_unc = self.choose_acceptance("batch", new_structure, new_structure_neighborlist, self.batched_error_constant, self.batched_temp, lmp_input, default_id)
            else: 
                keep_new, new_error = self.choose_acceptance("batch", new_structure, new_structure_neighborlist, self.batched_error_constant)

            if keep_new:
                self.success_step += 1
                new_error = round(new_error, 7)
                self.current_error = new_error
                error_list.append(f'{self.nsteps},{new_error}')
                moves += len(moved_atoms)

                if self.hybrid == True:
                    new_energy, max_unc = self.worker_task(new_structure, lmp_acc_in, default_id)
                    new_energy = round(new_energy, 5)
                    self.current_energy = new_energy
                    print(f'Step {self.nsteps}. Accepted, sum of residuals = {new_error}.  Energy = {new_energy/num_atoms} per atom')
                    with open('errors.txt', 'a') as out:
                        out.write(f"step # = {self.nsteps}, error = {new_error}, moved = {len(moved_atoms)}, tot_moves = {moves}, moves attempted = {moves_attempted}, error_const = {self.batched_error_constant}\n")
                    with open('error_plotting.txt', 'a') as out:
                        out.write(f'{self.nsteps} {new_error} {new_energy/num_atoms} {self.batched_error_constant} {self.batched_temp} {max_unc}\n')
                else:
                    print(f'Step {self.nsteps}. Accepted, sum of residuals = {new_error}.')
                    with open('errors.txt', 'a') as out:
                        out.write(f"step # = {self.nsteps}, error = {new_error}, moved = {len(moved_atoms)}, tot_moves = {moves}, moves attempted = {moves_attempted}, error_const = {self.batched_error_constant}\n")
                    with open('error_plotting.txt', 'a') as out:
                        out.write(f'{self.nsteps} {new_error} {self.batched_error_constant}\n')

                current_structure = copy.deepcopy(new_structure)                    
                current_structure_neighborlist = current_structure.get_all_neighbors(r=10.0)
                current_structure.to(fmt='POSCAR', filename='output.vasp')
                frame_num = len(error_list)
                current_structure.write_xdatcar(frame_num)
    
                counter.append(self.nsteps)

            if self.hybrid == True:
                #quenching scheme                
                if self.batched_temp > final_temp:
                    self.batched_temp = self.batched_temp * self.q_temp ** ( (self.nsteps * num_processes) / 1000000 )
                self.batched_error_constant = self.batched_error_constant *  self.q_scatter ** ( ( (self.nsteps * num_processes) / 1000000 ) / 2)
                #use this line to change the number in a batch. need rounding to return an integer value
                # num_batch = round(num_batch * batch_decay)
            
        output_structure = current_structure
        output_structure.to(fmt='POSCAR', filename='output.vasp')
        print(f'counter = {counter}')
        print(f'moves ={moves}')
        return output_structure
   

        
        
        
