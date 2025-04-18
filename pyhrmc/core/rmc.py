# -*- coding: utf-8 -*-
import os
import numpy as np
from numpy import exp
from random import random
import multiprocessing
import copy
from itertools import combinations
import warnings
import psutil, os
import logging

import pyhrmc.transformers as transform_mod
import pyhrmc.validators as validator_mod
from pyhrmc.core.structure import Structure
from pyhrmc.core.interpolator import Interpolator


logging.basicConfig(level=logging.INFO)


class RMC:
    def __init__(
        self,
        experimental_G_csv,
        sigma,
        q_scatter,
        q_temp=None,
        init_temp=None,
        hybrid=True,
        lmp_exec=None,
        dump_freq=5000
    ):
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
        self.dump_freq = dump_freq

        if self.hybrid == True:
            from pyhrmc.core.energy import Energy
            self.Lammps_HRMC = Energy
            self.lmp_exec = lmp_exec

        if self.hybrid == True and self.q_temp == None and self.batched_temp == None:
            raise RuntimeError(
                "If running HRMC, please provide values for q_temp and init_temp.\n"
                "If running RMC, specify that hybrid == False\n"
            )

    def apply_oxi_state(self, valences, struct):
        for i in range(len(struct.sites)):
            # just choose first set of oxi state guesses in the returned tuple
            struct.sites[i].oxi_state = valences[0][struct.sites[i].specie.symbol]
        return

    # generate and check function to parallelize
    def make_candidate_structure(
        self, current_structure, validator_objects, transformer, task_id, lmp_input=None
    ):
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
                keep_new, new_error, max_unc = self.choose_acceptance(
                    "single",
                    new_structure,
                    neighborlist,
                    self.batched_error_constant,
                    self.batched_temp,
                    lmp_input,
                    task_id,
                )
            else:
                keep_new = self.choose_acceptance(
                    "single", 
                    new_structure, 
                    neighborlist, 
                    self.batched_error_constant
                )
            if keep_new:
                pass
            else:
                new_structure = False
        return new_structure

    def choose_acceptance(
        self,
        version,
        structure,
        neighborlist,
        error_constant,
        temp=None,
        lmp_input=None,
        task_id=None,
        validator_objects=None
    ):
        max_unc = 0
        if self.hybrid == True:
            # HRMC

            lammps_run = self.Lammps_HRMC(lmp_input, task_id, self.lmp_exec)
            new_energy, max_unc = lammps_run.lammps_energy(
                structure, self.nsteps, self.success_step, self.dump_freq
            )
        structure.load_experimental_from_file(self.experimental_G_csv)
        new_error, slope = structure.prdf_error(neighborlist)

        keep_new = True

        """
        If batched, then check the Validators again with the new structure
        """

        if version == "batch":
            is_valid = True  # true until proven otherwise
            for validator in validator_objects:
                if not validator.check_structure(structure):
                    is_valid = False
                    break  # one failure is enough to stop and reject the structure
            if not is_valid:
                keep_new = False
                if self.hybrid == True:
                    return keep_new, new_error, max_unc
                elif self.hybrid == False:
                    return keep_new, new_error


        """
        REORGANIZE TO BE BASED ON IF HBRID == TRUE
        """

        if new_error < self.current_error:
            if self.hybrid == True:
                if new_energy <= self.current_energy:
                    # accept new structure
                    pass
                elif new_energy > self.current_energy:
                    # HRMC, with energy term
                    new_h_error = (new_error / (error_constant**2)) + (
                        new_energy / ((1.987204259 * 10**-3) * temp)
                    )
                    old_h_error = (self.current_error / (error_constant**2)) + (
                        self.current_energy / ((1.987204259 * 10**-3) * temp)
                    )
                    h_prob = exp(old_h_error - new_h_error)
                    if random() < h_prob:
                        pass
                    else:
                        keep_new = False
                    # logging probability info for record keeping
                return keep_new, new_error, max_unc
            # RMC
            elif self.hybrid == False:
                # accept new structure
                return keep_new, new_error

        else:
            if self.hybrid == True:
                # HRMC, with energy term
                new_h_error = (new_error / (error_constant**2)) + (
                    new_energy / ((1.987204259 * 10**-3) * temp)
                )
                old_h_error = (self.current_error / (error_constant**2)) + (
                    self.current_energy / ((1.987204259 * 10**-3) * temp)
                )
                h_prob = exp(old_h_error - new_h_error)
                if random() < h_prob:
                    pass
                else:
                    keep_new = False
                # logging probability info for record keeping
                return keep_new, new_error, max_unc
            elif self.hybrid == False:
                probability = exp(
                    -1 * (abs(new_error - self.current_error)) / error_constant
                )
                if random() < probability:
                    pass
                else:
                    keep_new = False
                return keep_new, new_error

    def worker_task(self, structure, lmp_input, task_id):
        lammps_run = self.Lammps_HRMC(lmp_input, task_id, self.lmp_exec)
        energy, max_unc = lammps_run.lammps_energy(
            structure, self.nsteps, self.success_step, self.dump_freq
        )
        return energy, max_unc

    def store_symbol(self, structure):
        el_switch = []
        j = 0
        el_list = [structure.sites[0].specie.symbol]
        for i in range(len(structure.sites)):
            if structure.atomic_numbers[i] != structure.atomic_numbers[j]:
                el_switch.append(i)
                el_list.append(structure.sites[i].specie.symbol)
            j = i
        return el_switch, el_list
    

    def coord_status(self, struct, validator_objects):
        validator = next((validator for validator in validator_objects if validator.__class__.__name__ == "DistancesCoordination"), None)
        if validator != None:
            composition = struct.composition.elements
            el_list = []
            for el_type in composition:
                el_list.append(el_type.symbol)
            for idx, site in enumerate(struct):
                # print(idx)
                points,sliced_df, voro = validator.get_voro(idx, struct)
                element_list, *_ = validator.get_coordination(idx, voro, sliced_df, points, struct)
                site.cn = {}
                for el_nn in el_list:
                    site.cn[el_nn] = len([nn for nn in element_list if nn == el_nn])
        return struct

    def run_rmc(
        self,
        num_processes=int(),
        initial_structure="amorphous_al2o3.vasp",
        experimental_G_csv="al2o3_5nm_gr.txt",
        keV=200,
        prdf_args={"bin_size": 0.04},
        transformations: dict = {
            "AtomHop": {},  # consider adding a second, smaller step size "max_step": 0.2
        },
        validators={},
        qmin = float(), 
        charges=None,
        # if defining charges as user:
        # charges = {
        # "species1": charge,
        # "species2" = charge,
        # ...
        # }
        TCS = None,
        pdf_cutoff=1.6,
        gaussian_blur = 1,
        lmp_init="in_init.lmp",
        lmp_test="in.lmp",
        lmp_accept="in_accept.lmp",
        max_steps=1000000,
        spec_order=None,
    ):

        if os.path.exists("error_plotting.txt"):
            if self.hybrid == True:
                try:
                    with open("last_step.txt", "r") as file:
                        restart_success_step = int(file.read())
                        self.success_step = restart_success_step
                except:
                    pass

                with open("error_plotting.txt", "r") as file:
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
                with open("error_plotting.txt", "r") as file:
                    lines = file.readlines()
                    last_line = lines[-1].strip()
                    variables = last_line.split()
                    self.nsteps = int(variables[0])
                    self.batched_error_constant = float(variables[2])

            warnings.warn(
                "\nA previous run has been detected. The simulation will resume at:\n"
                f"{self.nsteps} steps\n"
                "pyHRMC will automatically append the output of this run to existing output files.\n"
                "If this is undesired, please remove the existing file and restart the simulation."
            )

            if os.path.exists("XDATCAR"):
                warnings.warn("XDATCAR is missing.")


        """
        STAGING
        """

        # energy weighting staging
        final_temp = 100

        # batch number staging
        num_batch = num_processes
        # decay_const = 0.1
        # batch_decay = exp(- decay_const * num_batch)

        """
        INITIAL STRUCTURE AND CALCULATION SETUP
        """

        if self.hybrid == True:
            # lammps input setup
            lmp_init_in = []
            with open(lmp_init, "r") as file:
                for line in file:
                    lmp_init_in.append(line)

            lmp_input = []
            with open(lmp_test, "r") as file:
                for line in file:
                    lmp_input.append(line)

            lmp_acc_in = []
            with open(lmp_accept, "r") as file:
                for line in file:
                    lmp_acc_in.append(line)

        # convert to our python class
        initial_structure = Structure.from_file(initial_structure)
        initial_structure.xyz_df = initial_structure.xyz()
        # create neighborlist
        initial_structure_neighborlist = initial_structure.get_all_neighbors(r=10.0)
        num_atoms = len(initial_structure.sites)

        # cache/store interpolated radii
        struct_consts = Interpolator(initial_structure)


        if charges is not None:
            for site in initial_structure.sites:
                site.oxi_state = charges[site.specie.symbol]
            database_radii = struct_consts.database_ionic_radii()
            interpolated_radii = struct_consts.interpolated_ionic_radii(
                charges, database_radii
                )
            setattr(initial_structure, "interpolated_radii", interpolated_radii)
            print(f"Interpolated radii: {interpolated_radii}")
        else:
            charges = {}
            for species in initial_structure.types_of_species:
                el = species.name
                charges[el] = 0
            for site in initial_structure.sites:
                site.oxi_state = 0          


        if TCS == None:
            TCS = {}
            el_tuple = initial_structure.symbol_set
            for el in el_tuple:
                el_charge = charges[el]
                TCS[el] = struct_consts.interpolated_TCS(el, el_charge, keV, qmin)
        # if TCS != None, use user inputs 
        setattr(initial_structure, 'TCSs', TCS)
        print(f"Scattering cross sections: {TCS}")


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

        initial_structure = self.coord_status(initial_structure, validator_objects)

        """
        INITIALIZE RMC LOOP
        """
        current_structure = copy.deepcopy(initial_structure)
        current_structure.spec_order = spec_order
        el_switch, el_list = self.store_symbol(current_structure)
        current_structure.el_switch = el_switch
        current_structure.el_list = el_list
        current_structure.pdf_cutoff = pdf_cutoff
        current_structure.gaussian_blur = gaussian_blur
        current_structure_neighborlist = initial_structure_neighborlist
        self.current_error, slope = current_structure.prdf_error(
            current_structure_neighborlist
        )
        current_e = round(self.current_error, 5)
        # additional steps to initialize HRMC
        if self.hybrid == True:
            default_id = 100
            initial_e, max_unc = self.worker_task(
                current_structure, lmp_init_in, default_id
            )
            self.current_energy = initial_e
            print(
                f"Step {self.nsteps}. Sum of residuals = {current_e}. Energy = {self.current_energy/num_atoms}."
            )
        else:
            print(f"Step {self.nsteps}. Sum of residuals = {current_e}.")

        counter = []
        error_list = []
        moves = 0
        moves_attempted = 0

        current_structure.plot_pdf(
            current_structure_neighborlist, experimental_G_csv, slope
        )

        """
        RUN RMC LOOP
        """
        print("Starting RMC loops...")
        while self.nsteps < max_steps:
            self.nsteps += 1
            if self.nsteps % self.dump_freq == 0:
                if os.path.exists("pdfs.png"):
                    os.remove("pdfs.png")
                current_structure.plot_pdf(
                    current_structure_neighborlist, experimental_G_csv, slope
                )

            transformer = transformation_objects[
                0
            ]  # just grab the first transformation

            trial_structure = copy.deepcopy(current_structure)
            new_structure = copy.deepcopy(current_structure)

            with multiprocessing.Pool(processes=num_processes) as pool:
                tasks = list(range(num_batch))
                if self.hybrid == True:
                    inputs = [
                        (
                            trial_structure,
                            validator_objects,
                            transformer,
                            task_id,
                            lmp_input,
                        )
                        for task_id in tasks
                    ]
                else:
                    inputs = [
                        (trial_structure, validator_objects, transformer, task_id)
                        for task_id in tasks
                    ]
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
                # this becomes:
                distance = np.linalg.norm(i[1].coords - j[1].coords)
                try:  # only do this if min distance constraint is in use
                    min_distance = validator_objects[1].min_distances[
                        f"{i[1].species_string}", f"{j[1].species_string}"
                    ]
                    if distance < min_distance:
                        try:
                            moved_atoms.remove(i)
                        except:
                            pass
                except:
                    pass
            if len(moved_atoms) == 0:
                continue

            move_indices = []
            for c_idx, c_site in moved_atoms:
                new_structure.sites[c_idx] = c_site
                move_indices.append(c_idx)
            # edit neighbor list here
            # BUG: switch to individuvudal neighbor method for speedup
            new_structure_neighborlist = new_structure.get_all_neighbors(r=10.0)
            # update xyz_df
            new_structure.xyz_df = new_structure.xyz()
            #create list of all moved atoms in this batch
            new_structure.move_indices = move_indices

            if self.hybrid == True:
                keep_new, new_error, max_unc = self.choose_acceptance(
                    "batch",
                    new_structure,
                    new_structure_neighborlist,
                    self.batched_error_constant,
                    self.batched_temp,
                    lmp_input,
                    default_id,
                    validator_objects
                )
            else:
                keep_new, new_error = self.choose_acceptance(
                    "batch",
                    new_structure,
                    new_structure_neighborlist,
                    self.batched_error_constant,
                    validator_objects = validator_objects
                )

            if keep_new:
                self.success_step += 1
                new_error = round(new_error, 7)
                self.current_error = new_error
                error_list.append(f"{self.nsteps},{new_error}")
                moves += len(moved_atoms)

                if self.hybrid == True:
                    new_energy, max_unc = self.worker_task(
                        new_structure, lmp_acc_in, default_id
                    )
                    new_energy = round(new_energy, 5)
                    self.current_energy = new_energy
                    print(
                        f"Step {self.nsteps}. Accepted, sum of residuals = {new_error}.  Energy = {new_energy/num_atoms} per atom"
                    )
                    with open("error_plotting.txt", "a") as out:
                        out.write(
                            f"{self.nsteps} {new_error} {new_energy/num_atoms} {self.batched_error_constant} {self.batched_temp} {max_unc}\n"
                        )
                else:
                    print(
                        f"Step {self.nsteps}. Accepted, sum of residuals = {new_error}."
                    )
                    with open("error_plotting.txt", "a") as out:
                        out.write(
                            f"{self.nsteps} {new_error} {self.batched_error_constant}\n"
                        )

                current_structure = copy.deepcopy(new_structure)
                current_structure_neighborlist = current_structure.get_all_neighbors(
                    r=10.0
                )
                current_structure.to(fmt="POSCAR", filename="output.vasp")
                if self.success_step % 10 == 0:
                    frame_num = len(error_list)
                    current_structure.write_xdatcar(frame_num)

                counter.append(self.nsteps)

            if self.hybrid == True:
                # quenching scheme
                if self.batched_temp > final_temp:
                    self.batched_temp = self.batched_temp * self.q_temp ** (
                        (self.nsteps * num_processes) / 1000000
                    )
                    self.batched_error_constant = (
                        self.batched_error_constant
                        * self.q_scatter ** (((self.nsteps * num_processes) / 1000000) / 2)
                    )

        output_structure = current_structure
        output_structure.to(fmt="POSCAR", filename="output.vasp")
        print(f"counter = {counter}")
        print(f"moves ={moves}")
        return output_structure
