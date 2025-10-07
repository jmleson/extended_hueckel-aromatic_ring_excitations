import unittest
from fractions import Fraction

import sympy as sp

from IrreducibleRepresentation import IrreducibleRepresentation
from PointGroup import PointGroup
from SALC import SALC
from SymmetryOperation import SymmetryOperation


class PointGroupTest(unittest.TestCase):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)  # sehr wichtig!
        self.p = PointGroup(n=6)
        # revert Point group to C2v for 1,3-Butadien:
        self.p.n = 4
        # set_up_symmetry_operations:
        self.p.operations = []
        o = SymmetryOperation(n=self.p.n, name="E", transform_p=lambda i: i, amount=1)
        self.p.operations.append(o)
        o = SymmetryOperation(n=self.p.n, name="C2", transform_p=lambda i: [-4, -3, -2, -1][i - 1], amount=1)
        self.p.operations.append(o)
        o = SymmetryOperation(n=self.p.n, name="σv(xz)", transform_p=lambda i: [-1, -2, -3, -4][i - 1], amount=1)
        self.p.operations.append(o)
        o = SymmetryOperation(n=self.p.n, name="σv(yz)", transform_p=lambda i: [4, 3, 2, 1][i - 1], amount=1)
        self.p.operations.append(o)

        #set_up_irreducible_representations:
        self.p.irreducible_representations = []
        symmetry_names = [op.name for op in self.p.operations]
        i = IrreducibleRepresentation({j: 1 for j in symmetry_names}, name="A1")
        self.p.irreducible_representations.append(i)
        A2_chars = [1, 1, -1, -1]
        i = IrreducibleRepresentation(dict(zip(symmetry_names, A2_chars)), name="A2")
        self.p.irreducible_representations.append(i)
        B1_chars = [1, -1, 1, -1]
        i = IrreducibleRepresentation(dict(zip(symmetry_names, B1_chars)), name="B1")
        self.p.irreducible_representations.append(i)
        B2_chars = [1, -1, -1, 1]
        i = IrreducibleRepresentation(dict(zip(symmetry_names, B2_chars)), name="B2")
        self.p.irreducible_representations.append(i)


    def test_reducible_representations(self):
        reducible_representation = self.p.get_reducible_representation_for_ring_p_orbitals()
        #expected: {'E': 4, 'C2': 0, 'σv(xz)': -4, 'σv(yz)': 0}
        assert len(reducible_representation.keys()) == 4
        assert reducible_representation["E"] == 4
        assert reducible_representation["σv(xz)"] == -4
        assert reducible_representation["C2"] == 0
        assert reducible_representation["σv(yz)"] == 0

    def test_irreducible_representations(self):
        reducible_representation = {'E': 4, 'C2': 0, 'σv(xz)': -4, 'σv(yz)': 0}
        irreducible_representation = self.p.decomposing_into_irreducible_representations(reducible_representation)
        # expected: {'A1': 0, 'A2': 2, 'B1': 0, 'B2': 2}
        assert len(irreducible_representation.keys()) == 4
        assert irreducible_representation["A1"] == 0
        assert irreducible_representation["A2"] == 2
        assert irreducible_representation["B1"] == 0
        assert irreducible_representation["B2"] == 2

    def test_project(self):
        irreducible_representation = {'A1': 0, 'A2': 2, 'B1': 0, 'B2': 2}
        SALCs = self.p.project(irreducible_representation, p_orbital_index=1)
        # expecting A2 and B2 SALC
        assert len(SALCs) == 2

        for s in SALCs:
            p_symbols = sp.symbols(f"p1:{s.n + 1}")
            eq = sum(coeff * p_symbols[idx]
                     for idx, coeff in enumerate(s.prefactors_of_AOs, 0))
            assert sp.simplify(eq - s.equation) == 0

            if s.irred == "A2":
                assert s.prefactors_of_AOs == [1/2, 0, 0, -1/2]
            elif s.irred == "B2":
                assert s.prefactors_of_AOs == [1/2, 0, 0, +1/2]
            else:
                raise Exception("not supposed to be here")

        # 2nd p orbital:
        SALCs = self.p.project(irreducible_representation, p_orbital_index=2)
        # expecting A2 and B2 SALC
        assert len(SALCs) == 2

        for s in SALCs:
            p_symbols = sp.symbols(f"p1:{s.n + 1}")
            eq = sum(coeff * p_symbols[idx]
                     for idx, coeff in enumerate(s.prefactors_of_AOs, 0))
            assert sp.simplify(eq - s.equation) == 0

            if s.irred == "A2":
                assert s.prefactors_of_AOs == [0, 1/2, -1/2, 0]
            elif s.irred == "B2":
                assert s.prefactors_of_AOs == [0, 1/2, 1/2, 0]
            else:
                raise Exception("not supposed to be here")


    def test_SALC_norm(self):
        p1, p4 = sp.symbols(f"p1 p4")
        s = SALC(n=4, irred = "A2", p_orbital_prefactors = {p1: Fraction(1, 2), p4: Fraction(-1, 2)} )
        assert s.normalized == False
        s.norm()
        assert s.normalized == True
        expected_equation = (1/sp.sqrt(2)) * p1 -  (1/sp.sqrt(2)) * p4
        assert sp.simplify(s.equation - expected_equation) == 0
        assert s.prefactors_of_AOs == [1/sp.sqrt(2), 0, 0, -1/sp.sqrt(2)]





if __name__ == "__main__":
    unittest.main()