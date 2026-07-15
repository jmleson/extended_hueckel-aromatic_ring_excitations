import unittest
from fractions import Fraction

import sympy as sp

from src.main.MoleculeRepresentation import MoleculeRepresentation
from src.PointGroups.C2v import C2v
from src.main.SALC import SALC
from src.main.molecule_orbital import molecule_orbital
from src.solving.get_molecular_orbitals_from_secular_equation import get_molecular_orbitals_from_secular_equation
from testing.round_and_collect import round_and_collect


class TestLinearButadiene(unittest.TestCase):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.p = MoleculeRepresentation(n=4, n_instead_of_c=[], bound_cl_to_c_positions=[])
        ### revert Point group to C2v for 1,3-Butadien:
        self.p.pointgroup = C2v(n=4)
        self.p.circular = False
        ###############################################
        p1, p2, p3, p4, alpha, beta = sp.symbols(f"p1 p2 p3 p4 alpha beta")
        self.phi_1 = SALC(n=4, irred="B2", p_orbital_prefactors={p1: 1 / sp.sqrt(2), p4: 1 / sp.sqrt(2)}, orbital_symbol="p")
        self.phi_2 = SALC(n=4, irred="B2", p_orbital_prefactors={p2: 1 / sp.sqrt(2), p3: 1 / sp.sqrt(2)}, orbital_symbol="p")
        self.phi_3 = SALC(n=4, irred="A2", p_orbital_prefactors={p1: 1 / sp.sqrt(2), p4: -1 / sp.sqrt(2)}, orbital_symbol="p")
        self.phi_4 = SALC(n=4, irred="A2", p_orbital_prefactors={p2: 1 / sp.sqrt(2), p3: -1 / sp.sqrt(2)}, orbital_symbol="p")


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
        s = SALC(n=4, irred = "A2", p_orbital_prefactors = {p1: Fraction(1, 2), p4: Fraction(-1, 2)}, orbital_symbol="p" )
        assert s.normalized == False
        s.norm()
        assert s.normalized == True
        expected_equation = (1/sp.sqrt(2)) * p1 -  (1/sp.sqrt(2)) * p4
        assert sp.simplify(s.equation - expected_equation) == 0
        assert s.prefactors_of_AOs == [1/sp.sqrt(2), 0, 0, -1/sp.sqrt(2)]



    def test_h_eff_of_SALCs(self):
        self.phi_1.norm()
        self.phi_2.norm()
        self.phi_3.norm()
        self.phi_4.norm()
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

        for i in [self.phi_1, self.phi_2, self.phi_3, self.phi_4]:
            before = i.equation
            before_factors = i.prefactors_of_AOs
            i.norm()
            assert i.normalized == True
            assert i.prefactors_of_AOs == before_factors
            assert i.equation == before


    def test_get_effective_hamilton_matrix(self):
        alpha, beta = sp.symbols(f"alpha beta")
        self.phi_1.norm()
        self.phi_2.norm()
        self.phi_3.norm()
        self.phi_4.norm()
        # B2:
        h_matrix = self.p.get_effective_hamilton_matrix(SALCs = [self.phi_1, self.phi_2], irred="B2")
        expected = sp.Matrix([
            [ alpha,  beta         ],
            [ beta,   alpha + beta ]
        ])
        assert h_matrix.shape == expected.shape
        assert h_matrix.equals(expected)
        molecule_orbitals = get_molecular_orbitals_from_secular_equation(h_matrix, info="tst", sorting_dict_values={alpha: 0, beta: -1})
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
        molecule_orbitals = get_molecular_orbitals_from_secular_equation(h_matrix, info="tst", sorting_dict_values={alpha: 0, beta: -1})
        e_1 = alpha - beta * ((1 + sp.sqrt(5)) / 2)
        e_2 = alpha - beta * ((1 - sp.sqrt(5)) / 2)
        assert sp.simplify(e_2 - molecule_orbitals[0].eigenvalue) == 0
        assert sp.simplify(e_1 - molecule_orbitals[1].eigenvalue) == 0


    def test_get_energy_levels(self):
        alpha, beta = sp.symbols(f"alpha beta")
        round_to = 3
        result = self.p.get_energy_levels()

        assert len(result) == 4

        for (index, symmetry, equation) in [
            (0, "B2", alpha + 1.618 * beta),
            (1, "A2", alpha + 0.618 * beta),
            (2, "B2", alpha - 0.618 * beta),
            (3, "A2", alpha - 1.618 * beta)
        ]:
            assert isinstance(result[index], molecule_orbital)
            assert result[index].symmetry == symmetry
            assert sp.simplify(round_and_collect(result[index].eigenvalue, [alpha, beta], round_to)
                               - (equation)
                               ) == 0











if __name__ == "__main__":
    unittest.main()