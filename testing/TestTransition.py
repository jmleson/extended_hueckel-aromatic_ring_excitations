import unittest
import sympy as sp

from MoleculeRepresentation import MoleculeRepresentation
from MoleculeState import MoleculeState
from Transition import Transition
from molecule_orbital import molecule_orbital


class TestMoleculeRepresentation(unittest.TestCase):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        p = MoleculeRepresentation(n=4)
        p_mo_orbitals = p.get_energy_levels()
        p.set_to_s_orbitals()
        s_mo_orbitals = p.get_energy_levels()
        self.m = MoleculeState(bonding_p=p_mo_orbitals, antibonding_s=s_mo_orbitals)

        alpha, alpha_s, beta, beta_s = sp.symbols("alpha alpha_s beta beta_s")

        # self.get_first_test_case()


    def get_first_test_case(self):
        # phi_2
        mo = molecule_orbital()
        mo.eigenvalue = self.alpha
        mo.eigenvector = sp.Matrix([[1], [0]])
        mo.salcs = [{"factor": 1, "salc": self.phi_3}]
        mo.occupation = 0
        mo.symmetry = "Eg"

        mo_s = molecule_orbital()
        mo_s.eigenvalue = self.alpha_s
        mo_s.eigenvector = sp.Matrix([[1], [0]])
        mo_s.salcs = [{"factor": 1, "salc": self.phi_s_3}]
        mo_s.occupation = 0
        mo_s.symmetry = "Eu"

        self.m = MoleculeState(bonding_p=mo, antibonding_s=mo_s)
        return



    def get_occupation_expectations(self):
        alpha, alpha_s, beta, beta_s = sp.symbols("alpha alpha_s beta beta_s")
        diff_2 = alpha_s - alpha
        diff_3 = diff_2
        diff_4 = alpha_s - 2 * beta_s - alpha + 2 * beta
        return [
            # single occupation -> orbital energy / eigenvalue:
            {"occ": (1,0,0,0), "exp": alpha + 2 * beta, "transitioning energy": [], "expected": []},
            {"occ": (0, 1, 0, 0), "exp": alpha, "transitioning energy": [diff_2], "expected": [(0,1,0,0)]},
            {"occ": (0, 0, 1, 0), "exp": alpha, "transitioning energy": [diff_2], "expected": [(0,0,1,0)]},
            {"occ": (0, 0, 0, 1), "exp": alpha - 2 * beta, "transitioning energy": [diff_4], "expected": [(0,0,0,1)]},
            # 2 electrons:
            {"occ": (1,1,0,0) , "exp": 2*alpha + 2 * beta, "transitioning energy": [diff_2], "expected": [(0,1,0,0)]},
            {"occ":  (1, 0, 0, 1), "exp": 2 * alpha, "transitioning energy": [diff_4], "expected": [(0,0,0,1)]},
            {"occ": (2,0,0,0), "exp": 2 * alpha + 4 * beta, "transitioning energy": [], "expected": []},
            {"occ": (0, 2, 0, 0), "exp": 2*alpha, "transitioning energy": [diff_2], "expected": [(0,1,0,0)]},
            {"occ": (0, 0, 2, 0), "exp": 2*alpha, "transitioning energy": [diff_2], "expected": [(0,0,1,0)]},
            {"occ": (0, 0, 0, 2), "exp": 2 * alpha - 4 * beta, "transitioning energy": [diff_4], "expected": [(0,0,0,1)]},
            # 3 electrons:
            {"occ":(1, 1, 1, 0) , "exp": 3 * alpha + 2 * beta, "transitioning energy": [diff_2, diff_3],
                                                                "expected": [(0,1,0,0),(0,0,1,0)]},
            {"occ": (1, 1, 0, 1), "exp": 3 * alpha, "transitioning energy": [diff_2, diff_4],
                                                    "expected": [(0,1,0,0),(0,0,0,1)]},
            {"occ": (0, 1, 1, 1), "exp": 3 * alpha - 2 * beta , "transitioning energy": [ diff_2, diff_3, diff_4],
                                                                "expected": [(0,1,0,0),(0,0,1,0),(0,0,0,1)]},
            # 4 electrons:
            {"occ": (1, 1, 1, 1), "exp": 4 * alpha, "transitioning energy": [ diff_2, diff_3, diff_4],
                                                                "expected": [(0,1,0,0),(0,0,1,0),(0,0,0,1)]},
            {"occ": (2, 1, 1, 0), "exp": 4 * alpha + 4 * beta, "transitioning energy": [diff_2, diff_3],
                                                                "expected": [(0,1,0,0),(0,0,1,0)]},
            {"occ": (1, 2, 1, 0), "exp": 4 * alpha + 2 * beta, "transitioning energy": [diff_2, diff_3],
                                                                "expected": [(0,1,0,0),(0,0,1,0)]}
                        ]


    def test_calculate_energy_for_state_and_occupation(self):
        for occupation in self.get_occupation_expectations():
            assert self.m.calculate_energy_for_state_and_occupation(state=self.m.bonding_p,
                                                                occupation=occupation["occ"]) == occupation["exp"]

    def test_calculate_energy_before_transition(self):
        for occupation in self.get_occupation_expectations():
            self.m.set_occupation(p_occupation=occupation["occ"])
            assert self.m.calculate_energy_before_transition() == occupation["exp"]

    def test_get_changed_orbital_index(self):
        t = Transition(n=4, s_after_transition=(1,0,0,0), p_before_transition=(0,0,0,0))
        assert t.get_changed_orbital_index() == 0
        t = Transition(n=4, s_after_transition=(0, 1, 0, 0), p_before_transition=(0, 0, 0, 0))
        assert t.get_changed_orbital_index() == 1
        t = Transition(n=4, s_after_transition=(0, 0, 1, 0), p_before_transition=(0, 0, 0, 0))
        assert t.get_changed_orbital_index() == 2
        t = Transition(n=4, s_after_transition=(0, 0, 0, 1), p_before_transition=(0, 0, 0, 0))
        assert t.get_changed_orbital_index() == 3
        try:
            t = Transition(n=4, s_after_transition=(0, 1, 0, 1), p_before_transition=(0, 0, 0, 0))
            t.get_changed_orbital_index()
        except Exception:
            pass
        else:
            raise AssertionError("Expected an exception but none was raised")




    def test_get_thinkable_transitions(self):
        for occupation in self.get_occupation_expectations():
            print("OCC", occupation, flush=True)
            self.m.set_occupation(p_occupation=occupation["occ"])

            s_occupations = self.m.get_thinkable_transitions()
            assert len(s_occupations) == len(occupation["transitioning energy"])
            for s in range(len(s_occupations)):
                assert isinstance(s_occupations[s], Transition)
                assert s_occupations[s].n == 4
                assert s_occupations[s].p_before_transition == occupation["occ"]
                # print(occupation, s_occupations[s].s_occupation_after_transition,occupation["expected"])
                assert s_occupations[s].s_occupation_after_transition == occupation["expected"][s]
                print(occupation,"\n\t", s_occupations[s].get_transitioning_energy(),"==",occupation["transitioning energy"][s],"\n\t=",
                      s_occupations[s].orbital_to_excite_to.eigenvalue, "-", s_occupations[s].orbital_to_excite_of.eigenvalue,
                      )
                assert s_occupations[s].get_transitioning_energy() == occupation["transitioning energy"][s]
