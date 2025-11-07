
import math
import unittest
from fractions import Fraction
import sympy as sp

from MoleculeRepresentation import MoleculeRepresentation
from SALC import SALC, norm_and_group_SALCs
from is_multiple import is_multiple
from round_and_collect import round_and_collect


class TestPyrazine(unittest.TestCase):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)  # sehr wichtig!
        self.p = MoleculeRepresentation(n=6, bound_cl_to_c_positions=[], n_instead_of_c=[1,4])
        self.p.circular = True


    def test_reducible_representations(self):
        reducible_representation = self.p.get_reducible_representation_for_ring_p_orbitals()
        # expected: {'E': 6, 'C2(z)': 0, 'C2(y)': -2, 'C2(x)': 0, 'i': 0, 'σ(xy)': -6, 'σ(xz)': 0, 'σ(yz)': 2}
        print(reducible_representation)
        assert reducible_representation["E"] == 6
        assert reducible_representation["C2(y)"] == -2
        assert reducible_representation["σ(xy)"] == -6
        assert reducible_representation["σ(yz)"] == 2
        for other in reducible_representation.keys():
            if other not in ["E", "C2(y)", "σ(xy)", "σ(yz)"]:
                assert reducible_representation[other] == 0
        assert len(reducible_representation.items()) == 8
