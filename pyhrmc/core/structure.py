# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import itertools
from scipy.ndimage import gaussian_filter1d, zoom
from pymatgen.core import Structure
from pyhrmc.core.rdf import PDF

# from pymatgen.analysis.bond_valence import BVAnalyzer
import matplotlib.pyplot as plt
import re
from typing import Sequence
from pymatgen.core.sites import PeriodicSite, Site
from pymatgen.core.composition import Composition
import collections
import warnings

warnings.filterwarnings("ignore")

class Structure(Structure):
    """
    Adds a few extra methods to Structures that are assumed to be 2D slabs.

    The slab must be orthogonal to the z axis.
    """

    @property
    def thickness_z(self):
        max_z = max([coords[2] for coords in self.cart_coords])
        min_z = min([coords[2] for coords in self.cart_coords])
        thickness = max_z - min_z
        if self.lattice.abc[2] - thickness >= 5:
            return {"min_z": min_z, "max_z": max_z, "thickness": thickness}
        else:
            return None

    @property
    def slab_volume(self):

        # we need volume of our slab for the calculation of g(r)

        # BUG: assumes c vector along z axis
        if self.lattice.matrix[2][0] or self.lattice.matrix[2][1]:
            raise Exception(
                "This method requires the c lattice vector to be "
                "directly parallel to the z axis"
            )

        a = self.lattice.matrix[0]
        b = self.lattice.matrix[1]
        if self.thickness_z is None:
            c = self.lattice.matrix[2]
        else:
            c = [0, 0, self.thickness_z["thickness"]]
        volume = np.dot(np.cross(a, b), c)

        return volume

    @property
    def slab_density(self):
        return self.composition.weight / self.slab_volume

    @property
    def element_pairs(self):
        try:
            elements = self.el_list
        except:
            elements = self.symbol_set
        return list(itertools.combinations_with_replacement(elements, 2))

    def move_indices(self, move_indices):
        self.move_indices = move_indices

    def oxidation_state_list(structure):
        comp = structure.composition
        val = Composition(comp)
        valences = val.oxi_state_guesses()
        return valences

    def xyz(self):
        data = []
        for i, site in enumerate(self.sites):
            temp = [site.coords[0], site.coords[1], site.coords[2], site.species_string]
            data.append(temp)
        labels = ["x", "y", "z", "el"]

        df = pd.DataFrame(data, columns=labels)

        df_x = df.copy().sort_values("x")
        df_y = df.copy().sort_values("y")
        df_z = df.copy().sort_values("z")

        return {"df_x": df_x, "df_y": df_y, "df_z": df_z}

    def update_xyz(
        self, move_indices, new_position, move_vector, max_step, old_position
    ):
        """update position"""

        # old_position = copy.deepcopy(list(self.xyz_df['df_x'].loc[move_index][0:4]))

        for move_index in move_indices:
            print("before df update")
            print(self.xyz_df["df_x"].loc[move_index])

            for d in self.xyz_df.values():
                d.at[move_index, "x"] = new_position[0]
                d.at[move_index, "y"] = new_position[1]
                d.at[move_index, "z"] = new_position[2]

            print("after df update")
            print(self.xyz_df["df_x"].loc[move_index])

            """ sort df_x, df_y, df_z """

            df_x_original = self.xyz_df["df_x"].at[move_index, "x"]
            print(df_x_original)

            for i, df_name in enumerate(["df_x", "df_y", "df_z"]):
                df = self.xyz_df[df_name]
                axis = df_name[-1]
                df.sort_values(axis, inplace=True)

    # -------------------------------------------------------------------------
    # Methods/attributes for generating and caching pRDFs
    # -------------------------------------------------------------------------

    ##################
    # We need to now insert our pymatgen code that gets a fast RDF.
    ##################

    # Pymatgen-analysis-diffusion RDF utilities
    # https://github.com/materialsvirtuallab/pymatgen-analysis-diffusion/blob/master/pymatgen/analysis/diffusion/aimd/rdf.py
    #

    """
    Do we want Gaussian smearing?
    Do we divide by the number density rho to get g(r)?
    
    Some notes on the pymatgen.analysis.diffusion.aimd.rdf.RadialDistributionFunction
    The property .raw_rdf gives no smearing whereas .rdf does.
    Both .raw_rdf and .rdf  have the same integrated area.
    Both .raw_rdf and .rdf divide by rho (the number density).
    Both .raw_rdf and .rdf are single-counting rather than double counting.
    
    As implemented, we need to ensure the separation between slabs in z is
    greater than the RDF r_max.  Otherwise we accidentally count atoms in the
    vertical direction.
    
    We should therefore have a validator check that z_separation > r_max at the
    start of a calculation.
    
    How do we quickly generate a list of indices of just one atom type?
    Here, we call this atom_list_1, atom_list_2
    
    """

    @property
    def weightings(self):
        normalization = 0
        for el1, el2 in self.element_pairs:
            b1 = self.TCSs[el1]
            b2 = self.TCSs[el2]
            m1 = self.composition.get_atomic_fraction(el1)
            m2 = self.composition.get_atomic_fraction(el2)
            normalization += m1 * b1 * m2 * b2
        weightings = {}
        for pair in self.element_pairs:
            el1, el2 = pair
            b1 = self.TCSs[el1]
            b2 = self.TCSs[el2]
            m1 = self.composition.get_atomic_fraction(el1)
            m2 = self.composition.get_atomic_fraction(el2)
            numerator = m1 * b1 * m2 * b2
            weighting = numerator / normalization
            weightings[pair] = weighting

        return weightings

    def partial_rdfs(self, neighborlist):

        # OPTIMIZE: consider caching the creation/fitting of this class because
        # it never changes. (i.e move it to a separate method and use @cached_property)

        cutoff = 10
        prdf_maker = PDF(
            cutoff=cutoff, bin_size=0.04, el_switch=self.el_switch, el_list=self.el_list
        )
        prdf_maker.fit([self])

        # Above this line potentially cache
        # Below this line the function is slow.  Wasting time grabbing elements.

        # Grab the partial RDFs
        # bins, prdf_dict = self.compute_prdf(self)
        bins, prdf_dict = prdf_maker.compute_prdf(self, neighborlist)

        # replace the rdfs with smoothed rdfs
        # this part is optional / needs tweaking to optimize

        for pair, rdf in prdf_dict.items():
            rdf_normalized = rdf * self.volume / self.composition[pair[1]]
            rdf_smoothed = gaussian_filter1d(rdf_normalized, sigma= self.gaussian_blur)
            prdf_dict[pair] = rdf_smoothed

        return prdf_dict

    def full_pdf_g(self, neighborlist):

        rdfs_dict = self.partial_rdfs(neighborlist)
        weighting_dict = self.weightings

        # cutoff=10, bin_size=0.04 is hardcoded above --> array size will
        # always be 250.
        g = np.zeros(250)
        for pair in self.element_pairs:
            rdf = rdfs_dict[pair]
            weight = weighting_dict[pair]
            g += rdf * weight

        return g
    
    def define_u(self):
        bin_size= 0.04
        if self.thickness_z['thickness']:
            thickness = self.thickness_z['thickness']

        r_values = np.arange(0, 10, bin_size)

        u_values = np.where(
            r_values < thickness,
            r_values / (2 * thickness),
            1 - (thickness / (2 * r_values))
        )

        return u_values


    def full_pdf_G(self, neighborlist):

        # cutoff=10, bin_size=0.04 is hardcoded above --> array size will
        # always be 250. This info is also the "bins" variable in the
        # partial_rdfs method -- which we didn't save earlier
        r = np.arange(0, 10, 0.04)

        g = self.full_pdf_g(neighborlist)
        rho = self.num_sites / self.slab_volume  
        rho_rmc = self.num_sites / self.volume
        rho_correction = rho_rmc / rho

        if self.thickness_z:
            u = self.define_u()
            G = 4 * np.pi  * rho * r * (rho_correction * g + u - 1)
        else: 
            G = 4 * np.pi  * rho * r * (g - 1)
        return G

    # -------------------------------------------------------------------------
    # Methods for comparing to experiment
    # -------------------------------------------------------------------------

    # These values are empty until load_experimental_from_file is called
    r_experimental = None
    G_experimental = None

    def load_experimental_from_file(self, file_name):
        # NOTE: remove sep when using CSVs
        my_data = pd.read_csv(file_name, sep=",")

        # save as attribute so we don't need to reload later
        self.r_experimental = my_data.r.values
        self.G_experimental = my_data.gr.values

    def prdf_error(self, neighborlist):
        if self.r_experimental is None or self.G_experimental is None:
            raise Exception("Please call load_experimental_from_file first")

        exp = self.G_experimental
        calc = self.full_pdf_G(neighborlist)

        # we need to scale the calc pdf to be the same size (r) and same number
        # of bins (number of x values). Other option here to ensure the bin
        # size of the calc matches the experimental by default.
        # This code block is from...
        # https://stackoverflow.com/questions/55186867/
        current_size = len(calc)
        target_size = len(exp)
        zoom_factor = target_size / current_size
        calc_scaled = zoom(calc, zoom_factor)

        # we will zero the PDF at distances less than the user-specified cutoff, in Angstroms.

        zero_distance = self.pdf_cutoff
        PDF_maximum_radial_distance = 10

        zero_until = int(zero_distance / PDF_maximum_radial_distance * target_size)

        zero_array = np.zeros(zero_until)

        calc_scaled[:zero_until] = zero_array
        exp[:zero_until] = zero_array

        slope, error, _, _ = np.linalg.lstsq(
            calc_scaled.reshape(-1, 1), exp.reshape(-1, 1), rcond=None
        )

        error = error[0]
        slope = slope[0]

        return error, slope

    def plot_pdf(self, neighborlist, file_name, slope):
        r = self.r_experimental
        exp = self.G_experimental
        calc = self.full_pdf_G(neighborlist)

        current_size = len(calc)
        target_size = len(exp)
        zoom_factor = target_size / current_size
        calc_scaled = zoom(calc, zoom_factor) * slope

        zero_distance = self.pdf_cutoff
        PDF_maximum_radial_distance = 10
        zero_until = int(zero_distance / PDF_maximum_radial_distance * target_size)
        zero_array = np.zeros(zero_until)
        calc_scaled[:zero_until] = zero_array
        self.full_pdf_G_scaled = calc_scaled
        exp[:zero_until] = zero_array

        plt.plot(r, exp, label=f"{file_name}")
        plt.plot(r, calc_scaled, label="calculated")
        plt.legend()
        plt.savefig("pdfs.png")
        plt.show()
        plt.close()

        with open("pdf.txt", "w") as file:
            file.write(f"{' '.join(map(str, r))} \n{' '.join(map(str, calc_scaled))}\n")


        return

    def plot_error(self):
        try:
            f = np.loadtxt("error_plotting.txt", dtype=float)
            step = f[:, 0]
            error = f[:, 1]

            plt.plot(step, error, label="error")
            plt.legend()
            plt.savefig("errors.png")
            plt.show()
        except:
            pass

    def scaled_calc_pdf(self, slope):
        exp = self.G_experimental
        calc = self.full_pdf_G

        current_size = len(calc)
        target_size = len(exp)
        zoom_factor = target_size / current_size
        calc_scaled = zoom(calc, zoom_factor) * slope

        zero_distance = self.pdf_cutoff
        PDF_maximum_radial_distance = 10
        zero_until = int(zero_distance / PDF_maximum_radial_distance * target_size)
        zero_array = np.zeros(zero_until)
        calc_scaled[:zero_until] = zero_array
        self.full_pdf_G_scaled = calc_scaled
        exp[:zero_until] = zero_array

        return calc_scaled

    def write_xdatcar(self, frame_num):
        with open("XDATCAR", "a") as f:
            out = "from pyHRMC\n" "           1\n"
            out += (
                np.array2string(self.lattice.matrix, suppress_small=True)
                .replace("[", "")
                .replace("]", "")
                + "\n"
            )
            comp = self.composition.formula
            comp = re.split("(\d+)", comp)[:-1]
            el = " ".join(comp[::2]) + "\n"
            out += el
            out += " ".join(comp[1::2]) + "\n"
            out += f"Direct configuration = {frame_num}\n"

            for e in el.split():
                for i in range(len(self)):
                    site = self.sites[i]
                    if site.species_string == e:
                        out += (
                            np.array2string(site.frac_coords, suppress_small=True)
                            .replace("[", "")
                            .replace("]", "")
                            + "\n"
                        )
            f.write(out)

    # refactored version of pymatgen's get_all_neighbors method for speed up by vectorization
    # pymatgen is licensed under a MIT license, which can be found in LICENSE.md
    def get_all_neighbors(
        self,
        r: float,
        include_index: bool = False,
        include_image: bool = False,
        sites: Sequence[PeriodicSite] | None = None,
        numerical_tol: float = 1e-8,
    ):

        if sites is None:
            sites = self.sites
        center_indices, points_indices, images, distances = self.get_neighbor_list(
            r=r, sites=sites, numerical_tol=numerical_tol
        )
        if len(points_indices) < 1:
            return [[]] * len(sites)
        neighbor_dict: dict[int, list] = collections.defaultdict(list)
        atol = Site.position_atol
        all_sites = self.sites

        # create inividual masks for each check
        valid_neighbors_mask = distances > numerical_tol

        # Apply the mask and append the valid neighbors
        for cindex, pindex, d, is_valid in zip(
            center_indices, points_indices, distances, valid_neighbors_mask
        ):
            if is_valid:
                neighbor_dict[cindex].append({"nn_distance": d, "index": pindex})
            else:
                psite = all_sites[pindex]
                csite = sites[cindex]
                if (
                    # The or construct returns True immediately once one of the conditions are satisfied.
                    psite.species != csite.species
                    or (not np.allclose(psite.coords, csite.coords, atol=atol))
                    or (psite.properties != csite.properties)
                ):
                    neighbor_dict[cindex].append({"nn_distance": d, "index": pindex})

        neighbors = [neighbor_dict[i] for i in range(len(sites))]

        return neighbors
