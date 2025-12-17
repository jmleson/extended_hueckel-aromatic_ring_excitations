
import unittest
import sympy as sp

from MoleculeRepresentation import MoleculeRepresentation
from SALC import SALC, norm_and_group_SALCs
from solving.calculate_hueckel_secular_equation import calculate_hueckel_secular_equation
from is_multiple import is_multiple


class TestPyridine(unittest.TestCase):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)  # sehr wichtig!
        self.p = MoleculeRepresentation(n=6, bound_cl_to_c_positions=[], n_instead_of_c=[1])
        self.p.circular = True

        # SALCs, ! here not energetically sorted:
        p1, p2, p3, p4, p5, p6 = sp.symbols(f"p1 p2 p3 p4 p5 p6")
        self.phi_1 = SALC(n=7, irred="A2",
                          p_orbital_prefactors={p1: 0, p2: 1, p3: 0,
                                                p4: 0, p5: 0, p6: -1},
                          orbital_symbol="p")
        self.phi_2 = SALC(n=7, irred="A2",
                          p_orbital_prefactors={p1: 0, p2: 0, p3: 1,
                                                p4: 0, p5: -1, p6: 0},
                          orbital_symbol="p")
        self.phi_3 = SALC(n=7, irred="B2",
                          p_orbital_prefactors={p1: 1, p2: 0, p3: 0,
                                                p4: 0, p5: 0, p6: 0},
                          orbital_symbol="p")
        self.phi_4 = SALC(n=7, irred="B2",
                          p_orbital_prefactors={p1: 0, p2: 0, p3: 0,
                                                p4: 1, p5: 0, p6: 0},
                          orbital_symbol="p")
        self.phi_5 = SALC(n=7, irred="B2",
                          p_orbital_prefactors={p1: 0, p2: 1, p3: 0,
                                                p4: 0, p5: 0, p6: 1},
                          orbital_symbol="p")
        self.phi_6 = SALC(n=7, irred="B2",
                          p_orbital_prefactors={p1: 0, p2: 0, p3: 1,
                                                p4: 0, p5: 1, p6: 0},
                          orbital_symbol="p")



    def test_get_all_SALCs(self):
        self.phi_1.norm()
        self.phi_2.norm()
        self.phi_3.norm()
        self.phi_4.norm()
        self.phi_5.norm()
        self.phi_6.norm()
        # get_all_SALCs function:
        salcs = self.p.get_all_SALCs()
        assert len(salcs) == 6
        for s in salcs:
            s.norm()
            s.print()
            if s.irred == "A2":
                case1 = is_multiple(s.equation, self.phi_1.equation)#p2-p6
                case2 = is_multiple(s.equation, self.phi_2.equation)
                assert case1 or case2
            elif s.irred == "B2":
                case3 = is_multiple(s.equation, self.phi_3.equation)#p1
                case4 = is_multiple(s.equation, self.phi_4.equation)
                case5 = is_multiple(s.equation, self.phi_5.equation)
                case6 = is_multiple(s.equation, self.phi_6.equation)
                assert case3 or case4 or case5 or case6
            else:
                print(s.equation, s.irred)
                raise Exception("what is this?")


    def test_adjoint(self):
        alpha, beta, alpha_N, beta_N = sp.symbols(f"alpha beta alpha_N beta_N")
        # N:
        assert self.p.orbitals_adjoint(1,2) == beta_N
        assert self.p.orbitals_adjoint(1, 6) == beta_N
        assert self.p.orbitals_adjoint(1,1) == alpha_N
        for others in [3, 4, 5]:
            assert self.p.orbitals_adjoint(1, others) == 0

        # C2:
        assert self.p.orbitals_adjoint(2,2) == alpha
        assert self.p.orbitals_adjoint(2, 1) == beta_N
        assert self.p.orbitals_adjoint(2, 3) == beta
        for others in [4,5,6]:
            assert self.p.orbitals_adjoint(2, others) == 0

        #C4
        assert self.p.orbitals_adjoint(4, 4) == alpha
        for neighbor in [3, 5]:
            assert self.p.orbitals_adjoint(4, neighbor) == beta
        for others in [1, 2, 6]:
            assert self.p.orbitals_adjoint(4, others) == 0

    def test_get_effective_hamilton_matrix(self):
        alpha, beta, alpha_N, beta_N = sp.symbols(f"alpha beta alpha_N beta_N")
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
        for i in SALCs_by_irred["B2"]:
            i.print()
        h_matrix = self.p.get_effective_hamilton_matrix(SALCs=SALCs_by_irred["B2"], irred="B2")
        sp.pprint(h_matrix)
        expected = sp.Matrix([
            [alpha_N,               beta_N * sp.sqrt(2),    0,                  0],
            [beta_N * sp.sqrt(2),   alpha,                  beta,               0],
            [0,                     beta ,                  alpha,              sp.sqrt(2)*beta ],
            [0,                     0,                      sp.sqrt(2)*beta,    alpha ]
        ])
        assert h_matrix.shape == expected.shape
        assert h_matrix.equals(expected)

        calculate_hueckel_secular_equation(h_matrix, info="test", sorting_dict_values={alpha_N: 0, beta_N: -1, alpha: 0, beta: -1})


