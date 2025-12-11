import unittest
import sympy as sp

from MoleculeRepresentation import MoleculeRepresentation
from SALC import SALC
from TransitionIntegral import TransitionIntegral
from molecule_orbital import molecule_orbital


class TestTransitionIntegral(unittest.TestCase):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.p = MoleculeRepresentation(n=4, bound_cl_to_c_positions=[], n_instead_of_c=[])

        # P ORBITALS
        self.p_mo_orbitals = self.p.get_energy_levels()

        #S ORBITALS
        self.p.set_to_s_orbitals()
        self.s_mo_orbitals = self.p.get_energy_levels()

        self.p1, self.p2, self.p3, self.p4, self.s1, self.s2, self.s3, self.s4 = sp.symbols(f"p1 p2 p3 p4 s1 s2 s3 s4 ")
        self.alpha, self.beta, self.alpha_s, self.beta_s, self.delta = sp.symbols("alpha beta alpha_s beta_s delta")

        self.test_transition = {'s occupation after transition': (0, 1, 0, 0),
                                'Δ E by transition': -4*self.alpha + self.alpha_s - 4*self.beta}
        self.phi_1 = SALC(n=4, irred="A2u",
                          p_orbital_prefactors={self.p1: 1 / sp.sqrt(4), self.p2: 1 / sp.sqrt(4), self.p3: 1 / sp.sqrt(4),
                                                self.p4: 1 / sp.sqrt(4)},
                          orbital_symbol="p")
        self.phi_s_1 = SALC(n=4, irred="A1g",
                            p_orbital_prefactors={self.s1: 1 / sp.sqrt(4), self.s2: 1 / sp.sqrt(4), self.s3: 1 / sp.sqrt(4),
                                                  self.s4: 1 / sp.sqrt(4)},
                            orbital_symbol="s")
        self.phi_2 = SALC(n=4, irred="B2u",
                          p_orbital_prefactors={self.p1: 1 / sp.sqrt(4), self.p2: -1 / sp.sqrt(4), self.p3: 1 / sp.sqrt(4),
                                                self.p4: -1 / sp.sqrt(4)},
                          orbital_symbol="p")
        self.phi_s_2 = SALC(n=4, irred="B1g",
                            p_orbital_prefactors={self.s1: 1 / sp.sqrt(4), self.s2: -1 / sp.sqrt(4), self.s3: 1 / sp.sqrt(4),
                                                  self.s4: -1 / sp.sqrt(4)},
                            orbital_symbol="s")
        self.phi_3 = SALC(n=4, irred="Eg", p_orbital_prefactors={self.p1: 1 / sp.sqrt(2), self.p3: -1 / sp.sqrt(2)},
                          orbital_symbol="p")
        self.phi_s_3 = SALC(n=4, irred="Eg", p_orbital_prefactors={self.s1: 1 / sp.sqrt(2), self.s3: -1 / sp.sqrt(2)},
                            orbital_symbol="s")
        self.phi_4 = SALC(n=4, irred="Eg", p_orbital_prefactors={self.p2: 1 / sp.sqrt(2), self.p4: -1 / sp.sqrt(2)},
                          orbital_symbol="p")
        self.phi_s_4 = SALC(n=4, irred="Eg", p_orbital_prefactors={self.s2: 1 / sp.sqrt(2), self.s4: -1 / sp.sqrt(2)},
                            orbital_symbol="s")

    def get_first_test_case(self):
        mo = molecule_orbital()
        # mo.eigenvalue = self.alpha
        # mo.eigenvector = sp.Matrix([[1], [0]])
        mo.salcs = [{"factor": 1, "salc": self.phi_3}]
        # mo.occupation = 0
        # mo.symmetry = "Eg"

        mo_s = molecule_orbital()
        # mo_s.eigenvalue = self.alpha_s
        # mo_s.eigenvector = sp.Matrix([[1], [0]])
        mo_s.salcs = [{"factor": 1, "salc": self.phi_s_3}]
        # mo_s.occupation = 0
        # mo_s.symmetry = "Eu"

        return TransitionIntegral(bra=mo, ket=mo_s)

    def get_second_test_case(self):
        mo = molecule_orbital()
        # mo.eigenvalue = self.alpha
        # mo.eigenvector = sp.Matrix([[1]])
        mo.salcs = [{"factor": 3, "salc": self.phi_1}]
        # mo.occupation = 0
        # mo.symmetry = "A2u"

        mo_s = molecule_orbital()
        # mo_s.eigenvalue = self.alpha_s
        # mo_s.eigenvector = sp.Matrix([[1]])
        mo_s.salcs = [{"factor": 7, "salc": self.phi_s_1}]
        # mo_s.occupation = 0
        # mo_s.symmetry = "A1g"

        return TransitionIntegral(bra=mo, ket=mo_s)


    def test_get_multipliers_as_sum(self):
        t = self.get_first_test_case()

        bra_multipliers = t.get_multipliers_as_sum(t.bra)
        assert len(bra_multipliers) == 2
        assert bra_multipliers[0]["total_factor"] == 1 / sp.sqrt(2)
        assert bra_multipliers[1]["total_factor"] == -1 / sp.sqrt(2)
        for i in range(len(bra_multipliers)):
            assert bra_multipliers[i]["p_symbols"] == sp.Symbol(f'p{bra_multipliers[i]["p_number"]}')
        assert bra_multipliers[0]["p_symbols"] == self.p1
        assert bra_multipliers[1]["p_symbols"] == self.p3

        ket_multipliers = t.get_multipliers_as_sum(t.ket)
        assert len(ket_multipliers) == 2
        assert ket_multipliers[0]["total_factor"] == 1 / sp.sqrt(2)
        assert ket_multipliers[1]["total_factor"] == -1 / sp.sqrt(2)
        for i in range(len(ket_multipliers)):
            assert ket_multipliers[i]["p_symbols"] == sp.Symbol(f's{ket_multipliers[i]["p_number"]}')
        assert ket_multipliers[0]["p_symbols"] == self.s1
        assert ket_multipliers[1]["p_symbols"] == self.s3

        #
        t = self.get_second_test_case()

        bra_multipliers = t.get_multipliers_as_sum(t.bra)
        assert len(bra_multipliers) == 4
        for i in range(len(bra_multipliers)):
            assert bra_multipliers[i]["total_factor"] == 3 * 1 / sp.sqrt(4)
            assert bra_multipliers[i]["p_symbols"] == sp.Symbol(f'p{bra_multipliers[i]["p_number"]}')
        assert bra_multipliers[0]["p_symbols"] == self.p1
        assert bra_multipliers[1]["p_symbols"] == self.p2
        assert bra_multipliers[2]["p_symbols"] == self.p3
        assert bra_multipliers[3]["p_symbols"] == self.p4

        ket_multipliers = t.get_multipliers_as_sum(t.ket)
        assert len(ket_multipliers) == 4
        for i in range(len(ket_multipliers)):
            assert ket_multipliers[i]["total_factor"] == 7 * 1 / sp.sqrt(4)
            assert ket_multipliers[i]["p_symbols"] == sp.Symbol(f's{ket_multipliers[i]["p_number"]}')
        assert ket_multipliers[0]["p_symbols"] == self.s1
        assert ket_multipliers[1]["p_symbols"] == self.s2
        assert ket_multipliers[2]["p_symbols"] == self.s3
        assert ket_multipliers[3]["p_symbols"] == self.s4

    def test_multiply_out(self):
        t = self.get_first_test_case()
        result = t.multiply_out()
        assert result == self.delta

        t = self.get_second_test_case()
        result = t.multiply_out()
        assert result == 3*7*self.delta
