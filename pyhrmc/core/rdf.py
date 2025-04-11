# -*- coding: utf-8 -*-
"""
Edited matminer rdf.py for speed
Only the PartialRadialDistributionFunction class
This code is licensed under a BSD-style license that can be found in the LICENSE.md file

"""

import itertools
import math
from copy import copy

import numpy as np
from pymatgen.core.periodic_table import Element, Specie

from matminer.featurizers.base import BaseFeaturizer
import bisect


class PDF(BaseFeaturizer):
    """
    Compute the partial radial distribution functions (PRDF) of an xtal structure

    The PRDF of a crystal structure is the radial distribution function broken
    down for each pair of atom types.  The PRDF was proposed as a structural
    descriptor by [Schutt *et al.*]
    (https://journals.aps.org/prb/abstract/10.1103/PhysRevB.89.205118)

    Args:
        cutoff: (float) distance up to which to calculate the RDF.
        bin_size: (float) size of each bin of the (discrete) RDF.
        include_elems: (list of string), list of elements that must be included in PRDF
        exclude_elems: (list of string), list of elements that should not be included in PRDF

    Features:
        Each feature corresponds to the density of number of bonds
        for a certain pair of elements at a certain range of
        distances. For example, "Al-Al PRDF r=1.00-1.50" corresponds
        to the density of Al-Al bonds between 1 and 1.5 distance units
        By default, this featurizer generates RDFs for each pair
        of elements in the training set.

    """

    def __init__(
        self,
        cutoff=10.0,
        bin_size=0.1,
        el_switch=[],
        el_list=[],
        include_elems=(),
        exclude_elems=(),
    ):
        self.cutoff = cutoff
        self.bin_size = bin_size
        self.el_switch = el_switch
        self.el_list = el_list
        self.elements_ = None
        self.include_elems = list(
            include_elems
        )  # Makes sure the element lists are ordered
        self.exclude_elems = list(exclude_elems)

    def precheck(self, s):
        """
        Precheck the structure is ordered.
        Args:
            s: (pymatgen.Structure)
        Returns:
            (bool): True if passing precheck, false if failing
        """
        return s.is_ordered

    def fit(self, X, y=None):
        """Define the list of elements to be included in the PRDF. By default,
        the PRDF will include all of the elements in `X`

        Args:
            X: (numpy array nx1) structures used in the training set. Each entry
                must be Pymatgen Structure objects.
            y: *Not used*
            fit_kwargs: *not used*

        Returns:
            self
        """

        # Initialize list with included elements
        elements = {Element(e) for e in self.include_elems}

        # Get all of elements that appaer
        for strc in X:
            elements.update(
                [
                    e.element if isinstance(e, Specie) else e
                    for e in strc.composition.keys()
                ]
            )

        # Remove the elements excluded by the user
        elements.difference_update([Element(e) for e in self.exclude_elems])

        # Store the elements
        self.elements_ = [e.symbol for e in sorted(elements)]

        return self

    def featurize(self, s):
        """
        Get PRDF of the input structure.
        Args:
            s: Pymatgen Structure object.

        Returns:
            prdf, dist: (tuple of arrays) the first element is a
                    dictionary where keys are tuples of element
                    names and values are PRDFs.
        """

        if not s.is_ordered:
            raise ValueError("Disordered structure support not built yet")
        if self.elements_ is None:
            raise Exception("You must run 'fit' first!")

        dist_bins, prdf = self.compute_prdf(s)  # Assemble the PRDF for each pair

        # Convert the PRDF into a feature array
        zeros = np.zeros_like(dist_bins)  # Zeros if elements don't appear
        output = []
        for key in itertools.combinations_with_replacement(self.elements_, 2):
            output.append(prdf.get(key, zeros))

        # Stack them together
        return np.hstack(output)

    def compute_prdf(self, s, neighborlist):
        """Compute the PRDF for a structure

        Args:
            s: (Structure), structure to be evaluated
        Returns:
            dist_bins - float, start of each of the bins
            prdf - dict, where the keys is a pair of elements (strings),
                and the value is the radial distribution function for those paris of elements
        """

        # Get the composition of the array
        s = copy(s)
        s.remove_oxidation_states()
        composition = s.composition.fractional_composition.to_reduced_dict

        # Get the distances between all atoms, cutoff of sphere radius = 20.0
        neighbors_lst = neighborlist
        # neighbors_lst = s.get_all_neighbors(self.cutoff)
        # neighbors_lst = self.get_all_neighbors(r=20.0)
        # Sort neighbors by type
        distances_by_type = {}
        for p in itertools.product(composition.keys(), composition.keys()):
            distances_by_type[p] = []

        def get_symbol(s, site):
            idx = bisect.bisect_right(self.el_switch, site)
            return self.el_list[idx]

        for site, nlst in zip(
            enumerate(s.sites), neighbors_lst
        ):  # Each list is a list for each site
            my_elem = get_symbol(s, site[0])
            for neighbor in nlst:
                rij = neighbor["nn_distance"]
                n_elem = get_symbol(s, neighbor["index"])
                distances_by_type[(my_elem, n_elem)].append(rij)

        # Compute and normalize the prdfs
        prdf = {}
        dist_bins = self._make_bins()
        shell_volume = (
            4.0
            / 3.0
            * math.pi
            * (np.power(dist_bins[1:], 3) - np.power(dist_bins[:-1], 3))
        )
        for key, distances in distances_by_type.items():
            # Compute histogram of distances
            dist_hist, dist_bins = np.histogram(
                distances, bins=dist_bins, density=False
            )
            # Normalize
            n_alpha = composition[key[0]] * s.num_sites
            rdf = dist_hist / shell_volume / n_alpha

            prdf[key] = rdf

        return dist_bins[:-1], prdf

    def _make_bins(self):
        """Generate the edges of the bins for the PRDF

        Returns:
            [list of float], edges of the bins
        """
        return np.arange(0, self.cutoff + self.bin_size, self.bin_size)

    def feature_labels(self):
        if self.elements_ is None:
            raise Exception("You must run 'fit' first!")
        bin_edges = self._make_bins()
        labels = []
        for e1, e2 in itertools.combinations_with_replacement(self.elements_, 2):
            for r_start, r_end in zip(bin_edges, bin_edges[1:]):
                labels.append(f"{e1}-{e2} PRDF r={r_start:.2f}-{r_end:.2f}")
        return labels

    def citations(self):
        return [
            "@article{Schutt2014,"
            'author = {Sch{"{u}}tt, K. T. and Glawe, H. and Brockherde, F. '
            'and Sanna, A. and M{"{u}}ller, K. R. and Gross, E. K. U.},'
            "doi = {10.1103/PhysRevB.89.205118},"
            "journal = {Physical Review B},"
            "month = {may},number = {20},pages = {205118},"
            "title = {{How to represent crystal structures for machine learning:"
            " Towards fast prediction of electronic properties}},"
            "url = {http://link.aps.org/doi/10.1103/PhysRevB.89.205118},"
            "volume = {89},"
            "year = {2014}}"
        ]

    def implementors(self):
        return ["Logan Ward", "Saurabh Bajaj"]
