import unittest
from fractions import Fraction
import sympy as sp

from MoleculeRepresentation import MoleculeRepresentation
from SALC import SALC
from is_multiple import is_multiple
from molecule_orbital import molecule_orbital
from testing.round_and_collect import round_and_collect


class TestCyclobutadiene(unittest.TestCase):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)  # sehr wichtig!
        self.p = MoleculeRepresentation(n=4, bound_cl_to_c_positions=[], n_instead_of_c=[])
        self.p.circular = True

        ######
        p1, p2, p3, p4, s1, s2, s3, s4, alpha, beta = sp.symbols(f"p1 p2 p3 p4 s1 s2 s3 s4 alpha beta")
        self.phi_1 = SALC(n=4, irred="A2u",
                          p_orbital_prefactors={p1: 1 / sp.sqrt(4), p2: 1 / sp.sqrt(4), p3: 1 / sp.sqrt(4), p4: 1 / sp.sqrt(4)},
                          orbital_symbol="p")
        self.phi_s_1 = SALC(n=4, irred="A1g",
                          p_orbital_prefactors={s1: 1 / sp.sqrt(4), s2: 1 / sp.sqrt(4), s3: 1 / sp.sqrt(4),
                                                s4: 1 / sp.sqrt(4)},
                          orbital_symbol="s")
        self.phi_2 = SALC(n=4, irred="B2u",
                          p_orbital_prefactors={p1: 1 / sp.sqrt(4), p2: -1 / sp.sqrt(4), p3: 1 / sp.sqrt(4), p4: -1 / sp.sqrt(4)},
                          orbital_symbol="p")
        self.phi_s_2 = SALC(n=4, irred="B1g",
                          p_orbital_prefactors={s1: 1 / sp.sqrt(4), s2: -1 / sp.sqrt(4), s3: 1 / sp.sqrt(4),
                                                s4: -1 / sp.sqrt(4)},
                          orbital_symbol="s")
        self.phi_3 = SALC(n=4, irred="Eg", p_orbital_prefactors={p1: 1 / sp.sqrt(2), p3: -1 / sp.sqrt(2)},
                          orbital_symbol="p")
        self.phi_s_3 = SALC(n=4, irred="Eg", p_orbital_prefactors={s1: 1 / sp.sqrt(2), s3: -1 / sp.sqrt(2)},
                          orbital_symbol="s")
        self.phi_4 = SALC(n=4, irred="Eg", p_orbital_prefactors={p2: 1 / sp.sqrt(2), p4: -1 / sp.sqrt(2)},
                          orbital_symbol="p")
        self.phi_s_4 = SALC(n=4, irred="Eg", p_orbital_prefactors={s2: 1 / sp.sqrt(2), s4: -1 / sp.sqrt(2)},
                          orbital_symbol="s")

    def test_reducible_representations(self):
        reducible_representation = self.p.get_reducible_representation_for_ring_p_orbitals()
        # expected: {'E': 4, 'C4(z)': 0, 'C2': 0, "C'2": -2, "2C''2": 0, 'i': 0, 'S4': 0, 'σh': -4, 'σv': 2, 'σd': 0}
        assert len(reducible_representation.keys()) == 10
        assert reducible_representation["E"] == 4
        assert reducible_representation['C4(z)'] == 0
        assert reducible_representation['C2'] == 0
        assert reducible_representation["C'2"] == -2
        assert reducible_representation["C''2"] == 0
        assert reducible_representation['i'] == 0
        assert reducible_representation['S4'] == 0
        assert reducible_representation['σh'] == -4
        assert reducible_representation['σv'] == 2
        assert reducible_representation['σd'] == 0

    def test_irreducible_representations(self):
        reducible_representation = {'E': 4, 'C4(z)': 0, 'C2': 0, "C'2": -2, "2C''2": 0, 'i': 0, 'S4': 0, 'σh': -4, 'σv': 2, 'σd': 0}
        irreducible_representation = self.p.decomposing_into_irreducible_representations(reducible_representation)
        # expected: {'A1g': 0, 'A2g': 0, 'B1g': 0, 'B2g': 0, 'Eg': 1, 'A1u': 0, 'A2u': 1, 'B1u': 0, 'B2u': 1, 'Eu': 0}
        assert len(irreducible_representation.keys()) == 10
        assert irreducible_representation["Eg"] == 1
        assert irreducible_representation["A2u"] == 1
        assert irreducible_representation["B2u"] == 1
        #all others 0:
        total_other = sum(
            v for k, v in irreducible_representation.items()
            if k not in {"Eg", "A2u", "B2u"}
        )
        assert total_other == 0

    def test_SALC_generation(self):
        irreducible_representation = {'A1g': 0, 'A2g': 0, 'B1g': 0, 'B2g': 0, 'Eg': 1, 'A1u': 0, 'A2u': 1, 'B1u': 0, 'B2u': 1, 'Eu': 0}
        # function of expected number of SALCs per irred:
        expectation = self.p.get_expected_number_of_SALC_per_irreducible_representation(irreducible_representation)
        assert expectation["Eg"] == 2
        assert expectation["A2u"] == 1
        assert expectation["B2u"] == 1
        # all others 0:
        assert sum(v for k, v in expectation.items() if k not in {"Eg", "A2u", "B2u"}  ) == 0

        # project function:
        SALCs = self.p.project(irreducible_representation, p_orbital_index=1)
        for s in SALCs:
            p_symbols = sp.symbols(f"p1:{s.n + 1}")
            eq = sum(coeff * p_symbols[idx]
                     for idx, coeff in enumerate(s.prefactors_of_AOs, 0))
            assert sp.simplify(eq - s.equation) == 0

            if s.irred == "Eg":
                assert s.prefactors_of_AOs == [1/4, 0, -1/4, 0]
            elif s.irred == "B2u":
                assert s.prefactors_of_AOs == [1/4, -1/4, 1/4, -1/4]
            elif s.irred == "A2u":
                assert s.prefactors_of_AOs == [1/4, 1/4, 1/4, 1/4]
            else:
                raise Exception("not supposed to be here")

        # get_all_SALCs function:
        salcs = self.p.get_all_SALCs()
        assert len(salcs) == 4
        for s in SALCs:
            if s.irred == "A2u":
                assert is_multiple(s.equation, self.phi_1.equation)
            if s.irred == "B2u":
                assert is_multiple(s.equation, self.phi_2.equation)
            if s.irred == "Eg":
                assert is_multiple(s.equation, self.phi_3.equation) or is_multiple(s.equation, self.phi_4.equation)


    def test_SALC_norm(self):
        p1, p2, p3, p4 = sp.symbols(f"p1 p2 p3 p4")
        s = SALC(n=4, irred = "A2", p_orbital_prefactors = {p1: Fraction(1, 4), p2: Fraction(1, 4), p3: Fraction(1, 4), p4: Fraction(1, 4)}, orbital_symbol="p" )
        assert s.normalized == False
        s.norm()
        assert s.normalized == True
        expected_equation = (1/sp.sqrt(4)) * p1 + (1/sp.sqrt(4)) * p2 + (1/sp.sqrt(4)) * p3 + (1/sp.sqrt(4)) * p4
        assert sp.simplify(s.equation - expected_equation) == 0
        assert s.prefactors_of_AOs == [1/sp.sqrt(4), 1/sp.sqrt(4), 1/sp.sqrt(4), 1/sp.sqrt(4)]

        for i in [self.phi_1, self.phi_2, self.phi_3, self.phi_4]:
            before = i.equation
            before_factors = i.prefactors_of_AOs
            i.norm()
            assert i.normalized == True
            assert i.prefactors_of_AOs == before_factors
            assert i.equation == before
#


    def test_h_eff_of_SALCs(self):
        self.phi_1.norm()
        self.phi_2.norm()
        self.phi_3.norm()
        self.phi_4.norm()
        alpha, beta = sp.symbols(f"alpha beta")

        # A2u:
        h = self.p.h_eff(salc_1=self.phi_1, salc_2=self.phi_1)
        assert h == alpha + 2 * beta

        # B2u:
        h = self.p.h_eff(salc_1=self.phi_2, salc_2=self.phi_2)
        assert h == alpha - 2 * beta

        # Eg:
        h = self.p.h_eff(salc_1=self.phi_3, salc_2=self.phi_3)
        assert h == alpha

        h = self.p.h_eff(salc_1=self.phi_3, salc_2=self.phi_4)
        assert h == 0

        h = self.p.h_eff(salc_1=self.phi_4, salc_2=self.phi_3)
        assert h == 0

        h = self.p.h_eff(salc_1=self.phi_4, salc_2=self.phi_4)
        assert h == alpha

    def test_get_energy_levels(self):
        alpha, beta = sp.symbols(f"alpha beta")
        round_to = 3
        result = self.p.get_energy_levels()
        assert len(result) == 4

        for (index, symmetry, equation) in [
                                                (0, "A2u", alpha + 2 * beta ),
                                                (1, "Eg", alpha),
                                                (2, "Eg", alpha),
                                                (3, "B2u", alpha - 2 * beta)
                                            ]:
            assert isinstance(result[index], molecule_orbital)
            assert result[index].symmetry == symmetry
            assert sp.simplify(round_and_collect(result[index].eigenvalue, [alpha,beta], round_to)
                               - ( equation)
                               ) == 0

############################ S ORBITALS ###################################
    def test_S_reducible_representations(self):
        self.p.set_to_s_orbitals()
        reducible_representation = self.p.get_reducible_representation_for_ring_p_orbitals()
        # expected: absolute of p version: {'E': 4, 'C4(z)': 0, 'C2': 0, "C'2": 2, "2C''2": 0, 'i': 0, 'S4': 0, 'σh': 4, 'σv': 2, 'σd': 0}
        assert len(reducible_representation.keys()) == 10
        assert reducible_representation["E"] == 4
        assert reducible_representation['C4(z)'] == 0
        assert reducible_representation['C2'] == 0
        assert reducible_representation["C'2"] == 2
        assert reducible_representation["C''2"] == 0
        assert reducible_representation['i'] == 0
        assert reducible_representation['S4'] == 0
        assert reducible_representation['σh'] == 4
        assert reducible_representation['σv'] == 2
        assert reducible_representation['σd'] == 0

    def test_irreducible_representations(self):
        self.p.set_to_s_orbitals()
        reducible_representation = {'E': 4, 'C4(z)': 0, 'C2': 0, "C'2": 2, "C''2": 0, 'i': 0, 'S4': 0, 'σh': 4, 'σv': 2, 'σd': 0}
        irreducible_representation = self.p.decomposing_into_irreducible_representations(reducible_representation)
        # expected: {'A1g': 1, 'A2g': 0, 'B1g': 1, 'B2g': 0, 'Eg': 0, 'A1u': 0, 'A2u': 0, 'B1u': 0, 'B2u': 0, 'Eu': 1}
        assert len(irreducible_representation.keys()) == 10
        assert irreducible_representation["Eu"] == 1
        assert irreducible_representation["A1g"] == 1
        assert irreducible_representation["B1g"] == 1
        #all others 0:
        total_other = sum(
            v for k, v in irreducible_representation.items()
            if k not in {"Eu", "A1g", "B1g"}
        )
        assert total_other == 0

    def test_SALC_generation(self):
        self.p.set_to_s_orbitals()
        irreducible_representation = {'A1g': 1, 'A2g': 0, 'B1g': 1, 'B2g': 0, 'Eg': 0, 'A1u': 0, 'A2u': 0, 'B1u': 0, 'B2u': 0, 'Eu': 1}
        # function of expected number of SALCs per irred:
        expectation = self.p.get_expected_number_of_SALC_per_irreducible_representation(irreducible_representation)
        assert expectation["Eu"] == 2
        assert expectation["A1g"] == 1
        assert expectation["B1g"] == 1
        # all others 0:
        assert sum(v for k, v in expectation.items() if k not in {"Eu", "A1g", "B1g"}) == 0

        # project function:
        SALCs = self.p.project(irreducible_representation, p_orbital_index=1)
        for s in SALCs:
            p_symbols = sp.symbols(f"s1:{s.n + 1}")
            eq = sum(coeff * p_symbols[idx]
                     for idx, coeff in enumerate(s.prefactors_of_AOs, 0))
            assert sp.simplify(eq - s.equation) == 0

            if s.irred == "Eu":
                assert s.prefactors_of_AOs == [1 / 4, 0, -1 / 4, 0]
            elif s.irred == "B1g":
                assert s.prefactors_of_AOs == [1 / 4, -1 / 4, 1 / 4, -1 / 4]
            elif s.irred == "A1g":
                assert s.prefactors_of_AOs == [1 / 4, 1 / 4, 1 / 4, 1 / 4]
            else:
                raise Exception("not supposed to be here")

        # get_all_SALCs function:
        salcs = self.p.get_all_SALCs()
        assert len(salcs) == 4
        for s in salcs:
            if s.irred == "A1g":
                assert is_multiple(s.equation, self.phi_s_1.equation)
            if s.irred == "B1g":
                assert is_multiple(s.equation, self.phi_s_2.equation)
            if s.irred == "Eu":
                assert is_multiple(s.equation, self.phi_s_3.equation) or is_multiple(s.equation, self.phi_s_4.equation)


    def test_h_eff_of_SALCs(self):
        self.p.set_to_s_orbitals()
        self.phi_s_1.norm()
        self.phi_s_2.norm()
        self.phi_s_3.norm()
        self.phi_s_4.norm()
        alpha_s, beta_s = sp.symbols(f"alpha_s beta_s")

        h = self.p.h_eff(salc_1=self.phi_s_1, salc_2=self.phi_s_1)
        assert h == alpha_s + 2 * beta_s

        h = self.p.h_eff(salc_1=self.phi_s_2, salc_2=self.phi_s_2)
        assert h == alpha_s - 2 * beta_s

        h = self.p.h_eff(salc_1=self.phi_s_3, salc_2=self.phi_s_3)
        assert h == alpha_s

        h = self.p.h_eff(salc_1=self.phi_s_3, salc_2=self.phi_s_4)
        assert h == 0

        h = self.p.h_eff(salc_1=self.phi_s_4, salc_2=self.phi_s_3)
        assert h == 0

        h = self.p.h_eff(salc_1=self.phi_s_4, salc_2=self.phi_s_4)
        assert h == alpha_s


if __name__ == "__main__":
    unittest.main()