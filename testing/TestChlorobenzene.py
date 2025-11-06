import math
import unittest
from fractions import Fraction
import sympy as sp

from MoleculeRepresentation import MoleculeRepresentation
from SALC import SALC, norm_and_group_SALCs
from is_multiple import is_multiple
from round_and_collect import round_and_collect


class TestChlorobenzene(unittest.TestCase):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)  # sehr wichtig!
        self.p = MoleculeRepresentation(n=6, bound_cl_to_c_positions=[1])
        self.p.circular = True


        # SALCs, ! here not energetically sorted:
        p1, p2, p3, p4, p5, p6, p7 = sp.symbols(f"p1 p2 p3 p4 p5 p6 p7")
        self.phi_1 = SALC(n=7, irred="A2",
                          p_orbital_prefactors={p1: 1, p2: 0, p3: 0,
                                                p4: 0, p5: 0, p6: 0, p7: 0},
                          orbital_symbol="p")
        self.phi_2 = SALC(n=7, irred="A2",
                          p_orbital_prefactors={p1: 0, p2: 0, p3: 0,
                                                p4: 0, p5: 0, p6: 1, p7: 0},
                          orbital_symbol="p")
        self.phi_3 = SALC(n=7, irred="A2",
                          p_orbital_prefactors={p1: 0, p2: 0, p3: 0,
                                                p4: 0, p5: 0, p6: 0, p7: 1},
                          orbital_symbol="p")
        self.phi_4 = SALC(n=7, irred="B2",
                          p_orbital_prefactors={p1: 0, p2: 1/2, p3: 0,
                                                p4: 0, p5: 0, p6:1/2, p7: 0},
                          orbital_symbol="p")
        self.phi_5 = SALC(n=7, irred="B2",
                          p_orbital_prefactors={p1: 0, p2: 1/2, p3: 0,
                                                p4: 0, p5: 0, p6: -1/2, p7: 0},
                          orbital_symbol="p")
        self.phi_6 = SALC(n=7, irred="B2",
                          p_orbital_prefactors={p1: 0, p2: 0, p3: 1/2,
                                                p4: 0, p5: 1/2, p6: 0, p7: 0},
                          orbital_symbol="p")
        self.phi_7 = SALC(n=7, irred="B2",
                          p_orbital_prefactors={p1: 0, p2: 0, p3: 1/2,
                                                p4: 0, p5: -1/2, p6: 0, p7: 0},
                          orbital_symbol="p")

    def test_reducible_representations(self):
        reducible_representation = self.p.get_reducible_representation_for_ring_p_orbitals()
        # expected: {'E': 7, 'C2': -3, 'σv(xz)': -7, 'σv(yz)': 3}
        print(reducible_representation)
        assert reducible_representation["E"] == 7
        assert reducible_representation["C2"] == -3
        assert reducible_representation["σv(xz)"] == -7
        assert reducible_representation["σv(yz)"] == 3
        assert len(reducible_representation.items()) == 4

    def test_irreducible_representations(self):
        reducible_representation = {'E': 7, 'C2': -3, 'σv(xz)': -7, 'σv(yz)': 3}
        irreducible_representation = self.p.decomposing_into_irreducible_representations(reducible_representation)
        # expected: {'A1': 0, 'A2': 2, 'B1': 0, 'B2': 5}
        assert len(irreducible_representation.keys()) == 4
        assert irreducible_representation["A2"] == 2
        assert irreducible_representation["B2"] == 5
        # all others 0:
        total_other = sum(
            v for k, v in irreducible_representation.items()
            if k not in {"A2", "B2"}
        )
        assert total_other == 0

    def test_SALC_generation(self):
        irreducible_representation = {'A1': 0, 'A2': 2, 'B1': 0, 'B2': 5}
        # function of expected number of SALCs per irred:
        expectation = self.p.get_expected_number_of_SALC_per_irreducible_representation(irreducible_representation)
        assert expectation["A1"] == 0
        assert expectation["A2"] == 2
        assert expectation["B1"] == 0
        assert expectation["B2"] == 5

        # project function:
        p_symbols = sp.symbols(f"p1:{7+1}")
        for C_in_z_axis in [1, 4, 7]:
            SALCs = self.p.project(irreducible_representation, p_orbital_index=C_in_z_axis)
            for s in SALCs:
                eq = sum(coeff * p_symbols[idx]
                         for idx, coeff in enumerate(s.prefactors_of_AOs, 0))
                assert sp.simplify(eq - s.equation) == 0

                if s.irred == "A2":
                    assert s.equation == 0
                if s.irred == "B2":
                    assert s.equation == p_symbols[C_in_z_axis-1]
                else:
                    raise Exception(f"not supposed to be here {s.irred}")

        SALCs = self.p.project(irreducible_representation, p_orbital_index=2)
        for s in SALCs:
            eq = sum(coeff * p_symbols[idx]
                     for idx, coeff in enumerate(s.prefactors_of_AOs, 0))
            assert sp.simplify(eq - s.equation) == 0
            # sp.pprint(s.equation)
            if s.irred == "A2":
                expected = [0, 1/2, 0, 0, 0, -1/2, 0]
                assert s.prefactors_of_AOs == expected
            elif s.irred == "B2":
                expected = [0, 1 / 2, 0, 0, 0, 1/2, 0]
                assert s.prefactors_of_AOs == expected
            else:
                raise Exception(f"not supposed to be here {s.irred}")

        SALCs = self.p.project(irreducible_representation, p_orbital_index=3)
        for s in SALCs:
            eq = sum(coeff * p_symbols[idx]
                     for idx, coeff in enumerate(s.prefactors_of_AOs, 0))
            assert sp.simplify(eq - s.equation) == 0
            # sp.pprint(s.equation)
            if s.irred == "A2":
                expected = [0, 0, 1/2, 0, -1/2, 0, 0]
                assert s.prefactors_of_AOs == expected
            elif s.irred == "B2":
                expected = [0, 0, 1/2,  0, 1/2, 0, 0]
                assert s.prefactors_of_AOs == expected
            else:
                raise Exception(f"not supposed to be here {s.irred}")

        # get_all_SALCs function:
        salcs = self.p.get_all_SALCs()
        assert len(salcs) == 7
        SALCs_by_irred = norm_and_group_SALCs(salcs)
        assert len(SALCs_by_irred["A2"]) == 2
        assert len(SALCs_by_irred["B2"]) == 5

    def test_h_eff_of_SALCs(self):
        self.phi_1.norm()
        self.phi_2.norm()
        self.phi_3.norm()
        self.phi_4.norm()
        self.phi_5.norm()
        self.phi_6.norm()
        self.phi_7.norm()
        alpha, beta, alpha_Cl, beta_Cl = sp.symbols(f"alpha beta alpha_Cl beta_Cl")
        round_to = 6

        # A2:
        h = self.p.h_eff(salc_1=self.phi_1, salc_2=self.phi_1)
        result = round_and_collect(h, [alpha, beta], round_to)
        assert sp.simplify(result - (alpha)) == 0
        h = self.p.h_eff(salc_1=self.phi_2, salc_2=self.phi_2)
        result = round_and_collect(h, [alpha, beta], round_to)
        assert sp.simplify(result - (alpha)) == 0
        h = self.p.h_eff(salc_1=self.phi_3, salc_2=self.phi_3)
        result = round_and_collect(h, [alpha, beta], round_to)
        assert sp.simplify(result - (alpha_Cl)) == 0
        # B2:
        h = self.p.h_eff(salc_1=self.phi_4, salc_2=self.phi_4)
        result = round_and_collect(h, [alpha, beta], round_to)
        assert sp.simplify(result - (alpha)) == 0
        h = self.p.h_eff(salc_1=self.phi_5, salc_2=self.phi_5)
        result = round_and_collect(h, [alpha, beta], round_to)
        assert sp.simplify(result - (alpha)) == 0
        h = self.p.h_eff(salc_1=self.phi_6, salc_2=self.phi_6)
        result = round_and_collect(h, [alpha, beta], round_to)
        assert sp.simplify(result - (alpha)) == 0
        h = self.p.h_eff(salc_1=self.phi_7, salc_2=self.phi_7)
        result = round_and_collect(h, [alpha, beta], round_to)
        assert sp.simplify(result - (alpha)) == 0

        #mixed terms:
        h = self.p.h_eff(salc_1=self.phi_1, salc_2=self.phi_3)
        result = round_and_collect(h, [alpha, beta], round_to)
        assert sp.simplify(result - (beta_Cl)) == 0


    def test_adjoint(self):
        alpha, beta, alpha_Cl, beta_Cl = sp.symbols(f"alpha beta alpha_Cl beta_Cl")
        # Cl:
        assert self.p.orbitals_adjoint(7,1) == beta_Cl
        assert self.p.orbitals_adjoint(7, 7) == alpha_Cl
        for i in range(1,6):
            assert self.p.orbitals_adjoint(i+1,7) == 0

        # C1:
        assert self.p.orbitals_adjoint(1, 1) == alpha
        assert self.p.orbitals_adjoint(1, 7) == beta_Cl
        for neighbor in [2, 6]:
            assert self.p.orbitals_adjoint(1,neighbor) == beta
        for others in [3,4,5]:
            assert self.p.orbitals_adjoint(1,others) == 0

        # C2:
        assert self.p.orbitals_adjoint(2,2) == alpha
        for neighbor in [1,3]:
            assert self.p.orbitals_adjoint(2, neighbor) == beta
        for others in [4,5,6,7]:
            assert self.p.orbitals_adjoint(2, others) == 0

    def test_get_effective_hamilton_matrix(self):
        alpha, beta, alpha_Cl, beta_Cl = sp.symbols(f"alpha beta alpha_Cl beta_Cl")

        SALCs = self.p.get_all_SALCs()
        SALCs_by_irred = norm_and_group_SALCs(SALCs)

        # A2:
        h_matrix = self.p.get_effective_hamilton_matrix(SALCs=SALCs_by_irred["A2"], irred="A2")
        expected = sp.Matrix([
            [alpha, beta],
            [beta, alpha]
        ])
        assert h_matrix.shape == expected.shape
        assert h_matrix.equals(expected)

        # B2:
        h_matrix = self.p.get_effective_hamilton_matrix(SALCs=SALCs_by_irred["B2"], irred="B2")
        # sp.pprint(h_matrix)
        expected = sp.Matrix([
            [alpha,             sp.sqrt(2) * beta,  0,                  0,                  beta_Cl],
            [sp.sqrt(2) * beta, alpha,              beta,               0,                  0],
            [0,                 beta,               alpha,              sp.sqrt(2) * beta,  0],
            [0,                 0,                  sp.sqrt(2) * beta,  alpha,              0],#
            [beta_Cl,              0,                  0,                  0,                  alpha_Cl]
        ])
        assert h_matrix.shape == expected.shape
        assert h_matrix.equals(expected)
