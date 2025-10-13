import math
import unittest
from fractions import Fraction
import sympy as sp

from MoleculeRepresentation import MoleculeRepresentation
from SALC import SALC
from is_multiple import is_multiple
from molecule_orbital import molecule_orbital
from round_and_collect import round_and_collect


class TestBenzene(unittest.TestCase):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)  # sehr wichtig!
        self.p = MoleculeRepresentation(n=6)
        self.p.circular = True

        p1, p2, p3, p4, p5, p6 = sp.symbols(f"p1 p2 p3 p4 p5 p6")
        self.phi_5 = SALC(n=6, irred="A2u",
                          p_orbital_prefactors={p1: 1 / 6, p2: 1 / 6, p3: 1 / 6,
                                                p4: 1 / 6, p5: 1 / 6, p6: 1 / 6},
                          orbital_symbol="p")
        self.phi_6 = SALC(n=6, irred="B2g",
                          p_orbital_prefactors={p1: 1 / 6, p2: -1 / 6, p3: 1 / 6,
                                                p4: -1 / 6, p5: 1 / 6, p6: -1 / 6},
                          orbital_symbol="p")
        self.phi_1 = SALC(n=6, irred="E1g",
                          p_orbital_prefactors={p1: 1 / 4, p2: +1 / 8, p3: -1 / 8,
                                                p4: -1 / 4, p5: -1 / 8, p6: +1 / 8},
                          orbital_symbol="p")
        self.phi_2 = SALC(n=6, irred="E1g",
                          p_orbital_prefactors={p1: 0, p2: +1 / 4, p3: +1 / 4,
                                                p4: 0, p5: -1 / 4, p6: -1 / 4},
                          orbital_symbol="p")
        self.phi_3 = SALC(n=6, irred="E2u",
                          p_orbital_prefactors={p1: 2 / 8, p2: -1 / 8, p3: -1 / 8, p4: 2 / 8, p5: -1 / 8, p6: -1 / 8},
                          orbital_symbol="p")
        self.phi_4 = SALC(n=6, irred="E2u",
                          p_orbital_prefactors={p1: 0, p2: -1 / 4, p3: 1 / 4, p4: 0, p5: -1 / 4, p6: 1 / 4},
                          orbital_symbol="p")



    def test_reducible_representations(self):
        reducible_representation = self.p.get_reducible_representation_for_ring_p_orbitals()
        print(reducible_representation)
        # expected: {'E': 6, 'C6': 0, 'C3': 0, 'C2': 0, "C'2": -2, "C''2": 0, 'i': 0, 'S3': 0, 'S6': 0, 'σh': -6, 'σd': 0, 'σv': 2}
        assert len(reducible_representation.keys()) == 12
        assert reducible_representation["E"] == 6
        assert reducible_representation["C\'2"] == -2
        assert reducible_representation["σh"] == -6
        assert reducible_representation["σv"] == 2
        # all others 0:
        total_other = sum(
            v for k, v in reducible_representation.items()
            if k not in {"E", "C\'2", "σh", "σv"}
        )
        assert total_other == 0

    def test_irreducible_representations(self):
        reducible_representation = {'E': 6, 'C6': 0, 'C3': 0, 'C2': 0, "C'2": -2, "C''2": 0, 'i': 0, 'S3': 0, 'S6': 0, 'σh': -6, 'σd': 0, 'σv': 2}
        irreducible_representation = self.p.decomposing_into_irreducible_representations(reducible_representation)
        # expected: {'A1g': 0, 'A2g': 0, 'B1g': 0, 'B2g': 1, 'E1g': 1, 'E2g': 0, 'A1u': 0, 'A2u': 1, 'B1u': 0, 'B2u': 0, 'E1u': 0, 'E2u': 1}
        assert len(irreducible_representation.keys()) == 12
        assert irreducible_representation["E1g"] == 1
        assert irreducible_representation["E2u"] == 1
        assert irreducible_representation["B2g"] == 1
        assert irreducible_representation["A2u"] == 1
        #all others 0:
        total_other = sum(
            v for k, v in irreducible_representation.items()
            if k not in {"E1g", "E2u", "B2g", "A2u"}
        )
        assert total_other == 0

    def test_SALC_generation(self):
        irreducible_representation = {'A1g': 0, 'A2g': 0, 'B1g': 0, 'B2g': 1, 'E1g': 1, 'E2g': 0, 'A1u': 0, 'A2u': 1, 'B1u': 0, 'B2u': 0, 'E1u': 0, 'E2u': 1}
        # function of expected number of SALCs per irred:
        expectation = self.p.get_expected_number_of_SALC_per_irreducible_representation(irreducible_representation)
        assert expectation["E1g"] == 2
        assert expectation["E2u"] == 2
        assert expectation["A2u"] == 1
        assert expectation["B2g"] == 1
        # all others 0:
        assert sum(v for k, v in expectation.items() if k not in {"E1g", "E2u", "B2g", "A2u"}  ) == 0




        # project function:
        SALCs = self.p.project(irreducible_representation, p_orbital_index=1)
        for s in SALCs:
            p_symbols = sp.symbols(f"p1:{s.n + 1}")
            eq = sum(coeff * p_symbols[idx]
                     for idx, coeff in enumerate(s.prefactors_of_AOs, 0))
            assert sp.simplify(eq - s.equation) == 0

            if s.irred == "E1g":
                expected = [1/6, +1/12, -1/12, -1/6, -1/12, +1/12]
                expected_expr = sum(expected[i] * p_symbols[i] for i in range(len(expected)))
                assert is_multiple(s.equation, expected_expr)
                s.norm()
                expected = self.phi_1.prefactors_of_AOs
                expected_expr = sum(expected[i] * p_symbols[i] for i in range(len(expected)))
                assert is_multiple(s.equation, expected_expr)
            elif s.irred == "B2g":
                expected = self.phi_6.prefactors_of_AOs
                actual = s.prefactors_of_AOs
                assert all(math.isclose(float(a), float(e), rel_tol=1e-9) for a, e in zip(actual, expected))
            elif s.irred == "A2u":
                expected = self.phi_5.prefactors_of_AOs
                actual = s.prefactors_of_AOs
                assert all(math.isclose(float(a), float(e), rel_tol=1e-9) for a, e in zip(actual, expected))
            elif s.irred=="E2u":
                expected = self.phi_3.prefactors_of_AOs
                expected_expr = sum(expected[i] * p_symbols[i] for i in range(len(expected)))
                sp.pprint(expected_expr)
                sp.pprint(s.equation)
                assert is_multiple(s.equation, expected_expr)
            else:
                raise Exception(f"not supposed to be here {s.irred}")

        # get_all_SALCs function:
        salcs = self.p.get_all_SALCs()
        assert len(salcs) == 6
        for s in salcs:
            if s.irred == "A2u":
                assert is_multiple(s.equation, self.phi_5.equation)
            elif s.irred == "B2g":
                assert is_multiple(s.equation, self.phi_6.equation)
            elif s.irred == "E1g":
                assert is_multiple(s.equation, self.phi_1.equation) or is_multiple(s.equation, self.phi_2.equation)
            elif s.irred == "E2u":
                assert is_multiple(s.equation, self.phi_3.equation) or is_multiple(s.equation, self.phi_4.equation)
            else:
                s.norm()
                print(s.equation, s.irred)
                raise Exception("what is this?")


    def test_h_eff_of_SALCs(self):
        self.phi_1.norm()
        self.phi_2.norm()
        self.phi_3.norm()
        self.phi_4.norm()
        self.phi_5.norm()
        self.phi_6.norm()
        alpha, beta = sp.symbols(f"alpha beta")
        round_to = 6

        # A2u:
        h = self.p.h_eff(salc_1=self.phi_5, salc_2=self.phi_5)
        assert sp.simplify(round_and_collect(h, [alpha,beta], round_to) - ( alpha + 2 * beta)) == 0
        # B2g:
        h = self.p.h_eff(salc_1=self.phi_6, salc_2=self.phi_6)
        assert sp.simplify(round_and_collect(h, [alpha, beta], round_to) - (alpha - 2 * beta)) == 0
        # E1g:
        h = self.p.h_eff(salc_1=self.phi_1, salc_2=self.phi_1)
        assert sp.simplify(round_and_collect(h, [alpha, beta], round_to) - (alpha + beta)) == 0
        h = self.p.h_eff(salc_1=self.phi_2, salc_2=self.phi_2)
        assert sp.simplify(round_and_collect(h, [alpha, beta], round_to) - (alpha + beta)) == 0
        # E2u:
        h = self.p.h_eff(salc_1=self.phi_3, salc_2=self.phi_3)
        assert sp.simplify(round_and_collect(h, [alpha, beta], round_to) - (alpha - beta)) == 0
        h = self.p.h_eff(salc_1=self.phi_4, salc_2=self.phi_4)
        assert sp.simplify(round_and_collect(h, [alpha, beta], round_to) - (alpha - beta)) == 0

