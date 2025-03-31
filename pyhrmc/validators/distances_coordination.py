# -*- coding: utf-8 -*-

# This code is adapted from pymatgen.analysis.local_env.py
# pymatgen is licensed under a MIT license, which can be found in LICENSE.md

from pyhrmc.validators import Validator
import numpy as np
import pandas as pd
import math
from pymatgen.analysis.molecule_structure_comparator import CovalentRadius
import warnings
from scipy.spatial import Voronoi
from collections import namedtuple


NNData = namedtuple("NNData", ["all_nninfo", "cn_weights", "cn_nninfo"])


class DistancesCoordination(Validator):
    def __init__(
        self,
        BulkCoordinationRange,
        SurfaceCoordinationRange = None,
        SurfaceDistance = None,
        MinDistances = None
    ):
        self.min_distances = MinDistances
        self.coordination = BulkCoordinationRange
        self.surface_coordination = SurfaceCoordinationRange
        self.surface_distance = SurfaceDistance
        self.NNData = NNData

        if self.surface_coordination is not None:
            if self.surface_distance is None:
                raise RuntimeError(
                            "If using surface coordination constraints, please select a value for the surface distance."
                        )
        if self.surface_distance is not None:
            if self.surface_coordination is None:
                raise RuntimeError(
                            "If using a surface distance, please select surface coordination constraints."
                        )

    @staticmethod
    def vol_tetra(vt1, vt2, vt3, vt4):
        """
        Calculate the volume of a tetrahedron, given the four vertices of vt1,
        vt2, vt3 and vt4.
        Args:
            vt1 (array-like): coordinates of vertex 1.
            vt2 (array-like): coordinates of vertex 2.
            vt3 (array-like): coordinates of vertex 3.
            vt4 (array-like): coordinates of vertex 4.
        Returns:
            (float): volume of the tetrahedron.
        """
        vol_tetra = np.abs(np.dot((vt1 - vt4), np.cross((vt2 - vt4), (vt3 - vt4)))) / 6
        return vol_tetra

    @staticmethod
    def solid_angle(center, coords):
        """
        Helper method to calculate the solid angle of a set of coords from the
        center.
        Args:
            center (3x1 array): Center to measure solid angle from.
            coords (Nx3 array): List of coords to determine solid angle.
        Returns:
            The solid angle.
        """

        # Compute the displacement from the center
        r = [np.subtract(c, center) for c in coords]

        # Compute the magnitude of each vector
        r_norm = [np.linalg.norm(i) for i in r]

        # Compute the solid angle for each tetrahedron that makes up the facet
        #  Following: https://en.wikipedia.org/wiki/Solid_angle#Tetrahedron
        angle = 0
        for i in range(1, len(r) - 1):
            j = i + 1
            tp = np.abs(np.dot(r[0], np.cross(r[i], r[j])))
            de = (
                r_norm[0] * r_norm[i] * r_norm[j]
                + r_norm[j] * np.dot(r[0], r[i])
                + r_norm[i] * np.dot(r[0], r[j])
                + r_norm[0] * np.dot(r[i], r[j])
            )
            if de == 0:
                my_angle = 0.5 * np.pi if tp > 0 else -0.5 * np.pi
            else:
                my_angle = np.arctan(tp / de)
            angle += (my_angle if my_angle > 0 else my_angle + np.pi) * 2

        return angle

    def get_site_statistics(self, voro, site_idx, struct, sliced_df):
        site_statistics_dict = {}
        all_vertices = voro.vertices
        center_coords = voro.points[site_idx]

        for nn, vind in voro.ridge_dict.items():
            # Get only those that include the site in question
            if site_idx in nn:
                other_site = nn[0] if nn[1] == site_idx else nn[1]

                # Get the solid angle of the face
                facets = [all_vertices[i] for i in vind]
                angle = self.solid_angle(center_coords, facets)

                # Compute the volume of associated with this face
                volume = 0
                # Notes from original code:
                # qvoronoi returns vertices in CCW order, so I can akak
                # the face up in to segments (0,1,2), (0,2,3), ... to compute
                # its area where each number is a vertex size
                for j, k in zip(vind[1:], vind[2:]):
                    volume += self.vol_tetra(
                        center_coords,
                        all_vertices[vind[0]],
                        all_vertices[j],
                        all_vertices[k],
                    )

                # Compute the distance of the site to the face
                true_other_site = sliced_df.iloc[other_site].name

                face_dist = (
                    np.linalg.norm(center_coords - struct.sites[true_other_site].coords)
                    / 2
                )

                # Compute the area of the face (knowing V=Ad/3)
                face_area = 3 * volume / face_dist

                # Compute the normal of the facet
                normal = np.subtract(
                    struct.sites[true_other_site].coords, center_coords
                )
                normal /= np.linalg.norm(normal)

                # Store by face index
                site_statistics_dict[other_site] = {
                    "site": struct.sites[true_other_site],
                    "normal": normal,
                    "solid_angle": angle,
                    "volume": volume,
                    "face_dist": face_dist,
                    "area": face_area,
                    "n_verts": len(vind),
                    "true_index": true_other_site,
                    "index": other_site,
                }

        return site_statistics_dict

    def _extract_nn_info(self, structure, nns):
        """Given Voronoi NNs, extract the NN info in the form needed by NearestNeighbors
        Args:
            structure (Structure): Structure being evaluated
            nns ([dicts]): Nearest neighbor information for a structure
        Returns:
            (list of tuples (Site, array, float)): See nn_info
        """

        # Extract the NN info
        siw = []
        max_weight = max(nn["solid_angle"] for nn in nns.values())
        for nstats in nns.values():
            site = nstats["site"]

            nn_info = {
                "site": site,
                "image": (0, 0, 0),
                "weight": nstats["solid_angle"] / max_weight,
                "true_site_index": nstats["true_index"],
                "site_index": nstats["index"],
            }

            siw.append(nn_info)
        return siw

    def get_default_radius(self, site):
        """
        An internal method to get a "default" covalent/element radius
        Args:
            site: (Site)
        Returns:
            Covalent radius of element on site, or Atomic radius if unavailable
        """
        try:
            return CovalentRadius.radius[site.specie.symbol]
        except Exception:
            return site.specie.atomic_radius

    def _get_radius(self, struct, move_index):
        """
        An internal method to get the expected radius for a site with
        oxidation state.
        Args:
            site: (Site)
        Returns:
            Oxidation-state dependent radius: ionic, covalent, or atomic.
            Returns 0 if no oxidation state or appropriate radius is found.
        """
        struct.sites[0]
        site = struct[move_index]
        el = site.species_string
        oxi = struct.sites[move_index].oxi_state

        if oxi == 0:
            return self.get_default_radius(site)
        else:
            for el_type in struct.interpolated_radii:
                if el_type[0] == el:
                    ionic_radius = el_type[1]

                    return ionic_radius
        return ionic_radius

    @staticmethod
    def transform_to_length(nndata, length):
        """
        Given NNData, transforms data to the specified fingerprint length
        Args:
            nndata: (NNData)
            length: (int) desired length of NNData
        """

        if length is None:
            return nndata

        if length:
            for cn in range(length):
                if cn not in nndata.cn_weights:
                    nndata.cn_weights[cn] = 0
                    nndata.cn_nninfo[cn] = []

        return nndata

    @staticmethod
    def _semicircle_integral(dist_bins, idx):
        """
        An internal method to get an integral between two bounds of a unit
        semicircle. Used in algorithm to determine bond probabilities.
        Args:
            dist_bins: (float) list of all possible bond weights
            idx: (float) index of starting bond weight
        Returns:
            (float) integral of portion of unit semicircle
        """

        r = 1

        x1 = dist_bins[idx]
        x2 = dist_bins[idx + 1]
        if x1 > 1:
            diff = x1 - 1
            x1 = x1 - diff
            x2 = x2 - diff
            if x2 < 0.0:
                x2 = 0.0


        if dist_bins[idx] == 1:
            area1 = 0.25 * math.pi * r**2
        else:
            area1 = 0.5 * (
                (x1 * math.sqrt(r**2 - x1**2))
                + (r**2 * math.atan(x1 / math.sqrt(r**2 - x1**2)))
            )

        area2 = 0.5 * (
            (x2 * math.sqrt(r**2 - x2**2))
            + (r**2 * math.atan(x2 / math.sqrt(r**2 - x2**2)))
        )

        return (area1 - area2) / (0.25 * math.pi * r**2)

    def check_distance(self, move_index, voro, sliced_df, points):
        """ Check that all bonds are longer than minimum distance"""
        # get absolute row position of move_index in sliced_df
        # this corresponds to absolute row position in the points array
        true_move_index = sliced_df.index.get_loc(move_index)
        neighbors = [x for x in voro.ridge_points if true_move_index in x]
        # ridge_points is a pair of atoms; we remove the pair that is not the center atom
        neighbors = [a if b == true_move_index else b for a, b in neighbors]

        # get distances to neighbors
        neighbor_points = np.array(points[neighbors], dtype=np.float64)
        center_points = np.array(points[true_move_index], dtype=np.float64)
        distances = np.linalg.norm(
            center_points - neighbor_points, axis=1
        )

        # compare distances to allowed distances based on element identity
        # return false if the distance is too short
        center_element = sliced_df.at[move_index, "el"]
        for distance, neighbor in zip(distances, neighbors):
            neighbor_element = sliced_df.iat[
                neighbor, sliced_df.columns.get_loc("el")
            ]
            if self.min_distances is not None:
                try:
                    allowed_distance = self.min_distances[
                        (center_element, neighbor_element)
                    ]
                except:
                    allowed_distance = self.min_distances[
                        (neighbor_element, center_element)
                    ]
                if distance < allowed_distance:
                    return False
        return True


    def get_coordination(self, move_index, voro, sliced_df, points, struct):

        neighbor_list = []
        element_list = []
        # get absolute row position of move_index in sliced_df
        # this corresponds to absolute row position in the points array
        true_move_index = sliced_df.index.get_loc(move_index)
        neighbors = [x for x in voro.ridge_points if true_move_index in x]
        # ridge_points is a pair of atoms; we remove the pair that is not the center atom
        neighbors = [a if b == true_move_index else b for a, b in neighbors]

        # get distances to neighbors
        neighbor_points = np.array(points[neighbors], dtype=np.float64)
        center_points = np.array(points[true_move_index], dtype=np.float64)
        distances = np.linalg.norm(
            center_points - neighbor_points, axis=1
        )
        center_element = sliced_df.at[move_index, "el"]

        """  Check that all coordination numbers fall within a range """

        nns = self.get_site_statistics(voro, true_move_index, struct, sliced_df)
        nn_info = self._extract_nn_info(struct, nns)

        # we are ignoring samples that have explicity porosity.
        # see https://github.com/materialsproject/pymatgen/blob/56b8c965ea6a70b4970df1b5b41396f9a1fd5f77/pymatgen/analysis/local_env.py#L3952

        # we are ignoring weighting of solid angle by electronegativity difference
        # see https://github.com/materialsproject/pymatgen/blob/56b8c965ea6a70b4970df1b5b41396f9a1fd5f77/pymatgen/analysis/local_env.py#L3973

        nn = sorted(nn_info, key=lambda x: x["weight"], reverse=True)

        """distance_cutoffs: ([float, float]) - if not None, penalizes neighbor
        distances greater than sum of covalent radii plus
        distance_cutoffs[0]. Distances greater than covalent radii sum
        plus distance_cutoffs[1] are enforced to have zero weight."""

        distance_cutoffs = (0.5, 1)

        if self.surface_distance is not None:
            bottom_surf = struct.thickness_z["min_z"] + self.surface_distance
            top_surf = struct.thickness_z["max_z"] - self.surface_distance

        # adjust solid angle weights based on distance
        # get radius for the moved atom, which is always [0] in move_indices
        r1 = self._get_radius(struct, move_index)
        for entry in nn:
            r2 = self._get_radius(struct, entry["true_site_index"])
            if r1 > 0 and r2 > 0:
                d = r1 + r2
            else:
                warnings.warn(
                    "CrystalNN: cannot locate an appropriate radius, "
                    "covalent or atomic radii will be used, this can lead "
                    "to non-optimal results."
                )
                d = self.get_default_radius(
                    struct[move_index]
                ) + self.get_default_radius(entry["site"])

            entry_coords = np.array(points[entry["site_index"]], dtype=np.float64)
            dist = np.linalg.norm(center_points - entry_coords)
            dist_weight: float = 0

            cutoff_low = d + distance_cutoffs[0]
            cutoff_high = d + distance_cutoffs[1]

            if dist <= cutoff_low:
                dist_weight = 1
            elif dist < cutoff_high:
                dist_weight = (
                    math.cos(
                        (dist - cutoff_low) / (cutoff_high - cutoff_low) * math.pi
                    )
                    + 1
                ) * 0.5
            entry["weight"] = entry["weight"] * dist_weight
            try:
                if center_points[2] < bottom_surf or center_points[2] > top_surf:
                    # corrective factor to account for the fact that weight is max 0.5
                    # if the center atom is at the surface
                    # this correction is "conservative" since atom may not
                    # be strictly at the surface
                    entry["weight"] = entry["weight"] * 1.2       
            except:
                pass             

        # sort nearest neighbors from highest to lowest weight
        nn = sorted(nn, key=lambda x: x["weight"], reverse=True)

        # setting maximum coordination number as 15
        length = 15
        nndata = ""
        nearby_atoms = []
        if nn[0]["weight"] == 0:
            nn = [x for x in nn if x["weight"] > 0]
            
        else:
            for entry in nn:
                entry["weight"] = round(entry["weight"], 3)

            # remove entries with no weight
            nn = [x for x in nn if x["weight"] > 0]

            # get the transition distances, i.e. all distinct weights
            dist_bins: list[float] = []
            for entry in nn:
                if not dist_bins or dist_bins[-1] != entry["weight"]:
                    dist_bins.append(entry["weight"])
            dist_bins.append(0)

            # main algorithm to determine fingerprint from bond weights
            cn_weights = {}  # CN -> score for that CN
            cn_nninfo = {}  # CN -> list of nearneighbor info for that CN
            for idx, val in enumerate(dist_bins):
                if val != 0:
                    nn_info = []
                    for entry in nn:
                        if entry["weight"] >= val:
                            nn_info.append(entry)
                    cn = len(nn_info)
                    cn_nninfo[cn] = nn_info
                    cn_weights[cn] = self._semicircle_integral(dist_bins, idx)

            # add zero coord
            cn0_weight = 1.0 - sum(cn_weights.values())
            if cn0_weight > 0:
                cn_nninfo[0] = []
                cn_weights[0] = cn0_weight

            nearby_atoms=[entry['true_site_index'] for entry in nn]

            nndata = self.transform_to_length(
                self.NNData(nn, cn_weights, cn_nninfo), length
            )

            max_key = max(nndata.cn_weights, key=lambda k: nndata.cn_weights[k])
            nn = nndata.cn_nninfo[max_key]
            for entry in nn:
                entry["weight"] = 1

        for n in nn:
            index = sliced_df.iloc[n["site_index"]].name
            neighbor_list.append(index)
            el = sliced_df.iloc[n["site_index"]]["el"]
            element_list.append(el)

        return element_list, center_element, neighbor_list, nearby_atoms

    def get_voro(self, move_index, struct):

        d = struct.xyz_df
        slice_distance = 9  # in Angstroms, should equal two coordnation spheres

        """Slice the structure around the atom in question"""
        x_pos = d["df_x"].at[move_index, "x"]
        lower_x = x_pos - slice_distance
        upper_x = x_pos + slice_distance

        y_pos = d["df_y"].at[move_index, "y"]
        lower_y = y_pos - slice_distance
        upper_y = y_pos + slice_distance

        z_pos = d["df_z"].at[move_index, "z"]
        lower_z = z_pos - slice_distance
        upper_z = z_pos + slice_distance

        sliced_df_x = d["df_x"].loc[
            (lower_x < d["df_x"]["x"]) & (d["df_x"]["x"] < upper_x)
        ]
        sliced_df_y = d["df_y"].loc[
            (lower_y < d["df_y"]["y"]) & (d["df_y"]["y"] < upper_y)
        ]
        sliced_df_z = d["df_z"].loc[
            (lower_z < d["df_z"]["z"]) & (d["df_z"]["z"] < upper_z)
        ]

        # if the slice extends beyond the simulation cell:
        """ wrap the slice around the simulation box """
        if lower_x < 0:
            low_df_x = d["df_x"].loc[d["df_x"]["x"] > struct.lattice.a + lower_x].copy()
            low_df_x["x"] = low_df_x["x"] - struct.lattice.a
            sliced_df_x = pd.concat(
                [sliced_df_x, low_df_x]
                )
        else:
            low_df_x = None
        if lower_y < 0:
            low_df_y = d["df_y"].loc[d["df_y"]["y"] > struct.lattice.b + lower_y].copy()
            low_df_y["y"] = low_df_y["y"] - struct.lattice.b
            sliced_df_y = pd.concat(
                [sliced_df_y, low_df_y]
                )
        else:
            low_df_y = None
        if lower_z < 0:
            low_df_z = d["df_z"].loc[d["df_z"]["z"] > struct.lattice.c + lower_z].copy()
            low_df_z["z"] = low_df_z["z"] - struct.lattice.c
            sliced_df_z = pd.concat(
                [sliced_df_z, low_df_z]
                )
        else:
            low_df_z = None
        if upper_x > struct.lattice.a:
            upper_df_x = d["df_x"].loc[d["df_x"]["x"] < upper_x - struct.lattice.a].copy()
            upper_df_x["x"] = upper_df_x["x"] + struct.lattice.a
            sliced_df_x = pd.concat(
                [sliced_df_x, upper_df_x]
                )
        else:
            upper_df_x = None
        if upper_y > struct.lattice.b:
            upper_df_y = d["df_y"].loc[d["df_y"]["y"] < upper_y - struct.lattice.b].copy()
            upper_df_y["y"] = upper_df_y["y"] + struct.lattice.b
            sliced_df_y = pd.concat(
                [sliced_df_y, upper_df_y]
                )
        else:
            upper_df_y = None
        if upper_z > struct.lattice.c:
            upper_df_z = d["df_z"].loc[d["df_z"]["z"] < upper_z - struct.lattice.c].copy()
            upper_df_z["z"] = upper_df_z["z"] + struct.lattice.c
            sliced_df_z = pd.concat(
                [sliced_df_z, upper_df_z]
                )
        else:
            upper_df_z = None


        # Get common items in sliced dataframes
        # Find common indices
        common_indices = sliced_df_x.index.intersection(sliced_df_y.index).intersection(sliced_df_z.index)

        sliced_df = sliced_df_x[
            sliced_df_x.isin(sliced_df_y.to_dict("list"))
            & sliced_df_x.isin(sliced_df_z.to_dict("list"))
        ].dropna()

        # List of dfs to check
        image_dfs = [low_df_x, low_df_y, low_df_z, upper_df_x, upper_df_y, upper_df_z]
        common_image_indices = []
        for df in image_dfs:
            if df is not None:
                common_image_indices.extend(common_indices.intersection(df.index)) 

        filtered_df = pd.DataFrame(index=common_image_indices, columns=['x', 'y', 'z', 'el'])
        for idx in common_image_indices:
            for col in ['x', 'y', 'z']:
                # Get the values for this index and column from each DataFrame
                values = {
                    'sliced_df_x': sliced_df_x.at[idx, col],
                    'sliced_df_y': sliced_df_y.at[idx, col],
                    'sliced_df_z': sliced_df_z.at[idx, col]
                }
                
                value_counts = pd.Series(values).value_counts()
                # Select the value that is different
                if len(value_counts) > 1:  
                    differing_value = value_counts.idxmin() 
                else:
                    # If all values are the same, just choose one
                    differing_value = values['sliced_df_x']
                
                # Assign the differing value to the result df
                filtered_df.at[idx, col] = differing_value
                filtered_df.at[idx, 'el'] = sliced_df_x.at[idx, 'el']

        sliced_df = pd.concat(
            [sliced_df, filtered_df]
            )

        """ Call Voronoi function """
        points = sliced_df[["x", "y", "z"]].to_numpy()
        voro = Voronoi(points)

        return points, sliced_df, voro


    def check_structure(self,struct):
        move_indices = struct.move_indices
        for move_index in move_indices:
            points, sliced_df, voro = self.get_voro(move_index, struct)
            """ Checking coordination number """
            # run pymatgen-modified coordination number function on moved atom
            try:
                result = self.check_distance(move_index, voro, sliced_df, points)
                if result == False:
                    return False
                element_list, el, neighbor_list, nearby_atoms = self.get_coordination(
                    move_index=move_index, voro=voro, sliced_df=sliced_df, points=points, struct=struct
                )
            except:
                return False  # when distances are too short

            # check if the atom is near a surface
            if self.surface_coordination is not None:
                if (
                    struct.sites[move_index].z
                    < struct.thickness_z["min_z"] + self.surface_distance
                    or struct.sites[move_index].z
                    > struct.thickness_z["max_z"] - self.surface_distance
                ):
                    cn_constraints = self.surface_coordination[el]
                    for el_nn in cn_constraints:
                        el_cn = len([nn for nn in element_list if nn == el_nn])
                        current_el_cn = struct.sites[move_index].cn[el_nn]
                        if (
                            self.surface_coordination[el][el_nn][0]
                            <= el_cn
                            <= self.surface_coordination[el][el_nn][1]
                        ):
                            struct.sites[move_index].cn[el_nn] = current_el_cn
                            pass

                        else:
                            if not (
                            self.surface_coordination[el][el_nn][0]
                            <= current_el_cn
                            <= self.surface_coordination[el][el_nn][1]
                        ):
                                struct.sites[move_index].cn[el_nn] = current_el_cn
                                pass # the coordination number was already outside of range before move
                            else:
                                return False  # when coordination number isn't right
                else:
                    cn_constraints = self.coordination[el]
                    for el_nn in cn_constraints:
                        el_cn = len([nn for nn in element_list if nn == el_nn])
                        current_el_cn = struct.sites[move_index].cn[el_nn]
                        if (
                            self.coordination[el][el_nn][0]
                            <= el_cn
                            <= self.coordination[el][el_nn][1]
                        ):
                            struct.sites[move_index].cn[el_nn] = el_cn
                            pass
                        else:
                            if not (
                                self.coordination[el][el_nn][0]
                                <= current_el_cn
                                <= self.coordination[el][el_nn][1]
                            ):
                                struct.sites[move_index].cn[el_nn] = el_cn
                                pass # the coordination number was already outside of range before move
                            else:
                                return False  # when coordination number isn't right
            else:
                cn_constraints = self.coordination[el]
                for el_nn in cn_constraints:
                    el_cn = len([nn for nn in element_list if nn == el_nn])
                    current_el_cn = struct.sites[move_index].cn[el_nn]
                    if (
                        self.coordination[el][el_nn][0]
                        <= el_cn
                        <= self.coordination[el][el_nn][1]
                    ):
                        struct.sites[move_index].cn[el_nn] = el_cn
                        pass
                    else:
                        if not (
                        self.coordination[el][el_nn][0]
                        <= current_el_cn
                        <= self.coordination[el][el_nn][1]
                    ):
                            struct.sites[move_index].cn[el_nn] = el_cn
                            pass # the coordination number was already outside of range before move
                        else:
                            return False  # when coordination number isn't right
 
            # run pymatgen-modified coordination number function on neighbor atoms
            for atom in nearby_atoms:
                try:
                    result = self.check_distance(move_index, voro, sliced_df, points)
                    if result == False:
                        return False
                    element_list, el, _, _  = self.get_coordination(
                        atom, voro, sliced_df, points, struct
                    )
                except:
                    return False  # when distances are too short

                # check if the atom is near a surface

                if self.surface_coordination is not None:
                    if (
                        struct.sites[atom].z
                        < struct.thickness_z["min_z"] + self.surface_distance
                        or struct.sites[atom].z
                        > struct.thickness_z["max_z"] - self.surface_distance
                    ):
                        cn_constraints = self.surface_coordination[el]
                        for el_nn in cn_constraints:
                            el_cn = len([nn for nn in element_list if nn == el_nn])
                            current_el_cn = struct.sites[atom].cn[el_nn]
                            if (
                                self.surface_coordination[el][el_nn][0]
                                <= el_cn
                                <= self.surface_coordination[el][el_nn][1]
                            ):
                                struct.sites[atom].cn[el_nn] = el_cn
                                pass
                            else:
                                if not (
                                self.surface_coordination[el][el_nn][0]
                                <= current_el_cn
                                <= self.surface_coordination[el][el_nn][1]
                            ):
                                    struct.sites[atom].cn[el_nn] = el_cn
                                    pass # the coordination number was already outside of range before move
                                else:
                                    return False  # when coordination number isn't right
                    else:
                        cn_constraints = self.coordination[el]
                        for el_nn in cn_constraints:
                            el_cn = len([nn for nn in element_list if nn == el_nn])
                            current_el_cn = struct.sites[atom].cn[el_nn]
                            if (
                                self.coordination[el][el_nn][0]
                                <= el_cn
                                <= self.coordination[el][el_nn][1]
                            ):
                                struct.sites[atom].cn[el_nn] = el_cn
                                pass
                            else:
                                if not (
                                    self.coordination[el][el_nn][0]
                                    <= current_el_cn
                                    <= self.coordination[el][el_nn][1]
                                ):
                                    struct.sites[atom].cn[el_nn] = el_cn
                                    pass # the coordination number was already outside of range before move
                                else:
                                    return False  # when coordination number isn't right
                else:
                    cn_constraints = self.coordination[el]
                    for el_nn in cn_constraints:
                        el_cn = len([nn for nn in element_list if nn == el_nn])
                        current_el_cn = struct.sites[atom].cn[el_nn]
                        if (
                            self.coordination[el][el_nn][0]
                            <= el_cn
                            <= self.coordination[el][el_nn][1]
                        ):
                            struct.sites[atom].cn[el_nn] = el_cn
                            pass
                        else:
                            if not (
                            self.coordination[el][el_nn][0]
                            <= current_el_cn
                            <= self.coordination[el][el_nn][1]
                        ):
                                struct.sites[atom].cn[el_nn] = el_cn
                                pass # the coordination number was already outside of range before move
                            else:
                                return False  # when coordination number isn't right

        return True
