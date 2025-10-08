import unittest
from fractions import Fraction

import sympy as sp

from IrreducibleRepresentation import IrreducibleRepresentation
from PointGroup import PointGroup
from SALC import SALC
from SymmetryOperation import SymmetryOperation
from tst.molecule_orbital import molecule_orbital
from tst.solve_saekular_equation import calculate
from round_and_collect import round_and_collect


class PointGroupTest(unittest.TestCase):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)  # sehr wichtig!
        self.p = PointGroup(n=6)
        # revert Point group to C2v for 1,3-Butadien:
        self.p.n = 4
        self.p.circular = False
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


        #####################################
        p1, p2, p3, p4, alpha, beta = sp.symbols(f"p1 p2 p3 p4 alpha beta")
        self.phi_1 = SALC(n=4, irred="B2", p_orbital_prefactors={p1: 1 / sp.sqrt(2), p4: 1 / sp.sqrt(2)})
        self.phi_2 = SALC(n=4, irred="B2", p_orbital_prefactors={p2: 1 / sp.sqrt(2), p3: 1 / sp.sqrt(2)})
        self.phi_3 = SALC(n=4, irred="A2", p_orbital_prefactors={p1: 1 / sp.sqrt(2), p4: -1 / sp.sqrt(2)})
        self.phi_4 = SALC(n=4, irred="A2", p_orbital_prefactors={p2: 1 / sp.sqrt(2), p3: -1 / sp.sqrt(2)})


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



    def test_h_eff_of_SALCs(self):
        alpha, beta = sp.symbols(f"alpha beta")

        h = self.p.h_eff(salc_1=self.phi_1, salc_2=self.phi_1)
        assert h == alpha

        h = self.p.h_eff(salc_1=self.phi_2, salc_2=self.phi_2)
        assert h == alpha + beta

        h = self.p.h_eff(salc_1=self.phi_1, salc_2=self.phi_2)
        assert h == beta

        h = self.p.h_eff(salc_1=self.phi_3, salc_2=self.phi_3)
        assert h == alpha

        h = self.p.h_eff(salc_1=self.phi_4, salc_2=self.phi_4)
        assert h == alpha - beta

        h = self.p.h_eff(salc_1=self.phi_3, salc_2=self.phi_4)
        assert h == beta


    def test_get_effective_hamilton_matrix(self):
        alpha, beta = sp.symbols(f"alpha beta")
        # B2:
        h_matrix = self.p.get_effective_hamilton_matrix(SALCs = [self.phi_1, self.phi_2], irred="B2")
        expected = sp.Matrix([
            [ alpha,  beta         ],
            [ beta,   alpha + beta ]
        ])
        assert h_matrix.shape == expected.shape
        assert h_matrix.equals(expected)
        molecule_orbitals = calculate(h_matrix, info="tst", sorting_dict_values={alpha: 0, beta: -1})
        e_1 = alpha + beta * ((1+sp.sqrt(5))/2)
        e_2 = alpha + beta * ((1-sp.sqrt(5))/2)
        assert sp.simplify(e_1 - molecule_orbitals[0].eigenvalue) == 0
        assert sp.simplify(e_2 - molecule_orbitals[1].eigenvalue) == 0

        # A2:
        h_matrix = self.p.get_effective_hamilton_matrix(SALCs=[self.phi_3, self.phi_4], irred="A2")
        expected = sp.Matrix([
            [alpha, beta],
            [beta, alpha - beta]
        ])
        assert h_matrix.shape == expected.shape
        assert h_matrix.equals(expected)
        molecule_orbitals = calculate(h_matrix, info="tst", sorting_dict_values={alpha: 0, beta: -1})
        e_1 = alpha - beta * ((1 + sp.sqrt(5)) / 2)
        e_2 = alpha - beta * ((1 - sp.sqrt(5)) / 2)
        assert sp.simplify(e_2 - molecule_orbitals[0].eigenvalue) == 0
        assert sp.simplify(e_1 - molecule_orbitals[1].eigenvalue) == 0


    def test_get_energy_levels(self):
        alpha, beta = sp.symbols(f"alpha beta")
        round_to = 3
        result = self.p.get_energy_levels()

        assert len(result) == 4

        # testing lowest energy level:
        assert isinstance(result[0], molecule_orbital)
        assert result[0].symmetry == "B2"
        assert sp.simplify(round_and_collect(result[0].eigenvalue, [alpha,beta], round_to)
                           - ( alpha + 1.618 * beta )
                           ) == 0
        # testing 2nd lowest energy level:
        assert isinstance(result[1], molecule_orbital)
        assert result[1].symmetry == "A2"
        assert sp.simplify(
                        round_and_collect(result[1].eigenvalue, [alpha,beta], round_to)
                        - (alpha + 0.618 * beta)
                        ) == 0
        # testing 2nd highest energy level:
        assert isinstance(result[2], molecule_orbital)
        assert result[2].symmetry == "B2"
        assert sp.simplify(
                        round_and_collect(result[2].eigenvalue, [alpha, beta], round_to)
                        - (alpha - 0.618 * beta)
                        ) == 0
        # testing highest energy level:
        assert isinstance(result[3], molecule_orbital)
        assert result[3].symmetry == "A2"
        assert sp.simplify(
                        round_and_collect(result[3].eigenvalue, [alpha, beta], round_to)
                        - (alpha - 1.618 * beta)
                        ) == 0









if __name__ == "__main__":
    unittest.main()