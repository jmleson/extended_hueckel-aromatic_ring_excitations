import unittest
import sympy as sp

from MoleculeState import MoleculeState
from Transition import Transition
from molecule_orbital import molecule_orbital


class TestMoleculeState(unittest.TestCase):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.m = MoleculeState(n=4, bound_cl_to_c_positions=[], n_instead_of_c=[])


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
            part_1, msg = self.m.calculate_energy_for_state_and_occupation(state=self.m.bonding_p,
                                                                occupation=occupation["occ"])
            assert part_1 == occupation["exp"]

    def test_calculate_energy_before_transition(self):
        for occupation in self.get_occupation_expectations():
            self.m.set_occupation(p_occupation=occupation["occ"])
            part, msg = self.m.calculate_energy_before_transition()
            assert part == occupation["exp"]

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
            self.m.set_occupation(p_occupation=occupation["occ"])

            s_occupations, msg = self.m.get_thinkable_transitions()
            assert len(s_occupations) == len(occupation["transitioning energy"])
            for s in range(len(s_occupations)):
                assert isinstance(s_occupations[s], Transition)
                assert s_occupations[s].n == 4
                assert s_occupations[s].p_before_transition == occupation["occ"]
                assert s_occupations[s].s_occupation_after_transition == occupation["expected"][s]
                assert s_occupations[s].get_transitioning_energy() == occupation["transitioning energy"][s]


    def test_get_p_after_transition(self):
        p_before = (3, 7, 4, 2, 3, 10)
        with self.assertRaises(Exception):
            self.m.set_occupation(p_occupation=p_before)

        p_before = (4, 2, -3, 10)
        with self.assertRaises(Exception):
            self.m.set_occupation(p_occupation=p_before)

        # Normale Funktionalität testen
        p_before = (2, 1, 1, 0)
        self.m.set_occupation(p_occupation=p_before)
        result = self.m.get_p_after_transition(s_after_transition=(1, 0, 0, 0))
        assert result == (1, 1, 1, 0)
        result = self.m.get_p_after_transition(s_after_transition=(0, 1, 0, 0))
        assert result == (2, 0, 1, 0)

    def test_symmetry_allowed_transition(self):
        t = Transition(n=4, s_after_transition=(1,0,0,0), p_before_transition=(0,0,0,0))
        mo_p = molecule_orbital()
        mo_s = molecule_orbital()


        for combi in [("Eg", "Eu", True), ("A2u", "A1g", True), ("B2u", "B1g", True),
                      ("Eg", "A1g", False), ("Eg", "B2u", False), ("Eg", "B1g", False),
                      ("A1g", "A2u", True), ("A1g", "B2u", False), ("A1g", "A1g", False),
                      ]:
            mo_p.symmetry = combi[0]
            mo_s.symmetry = combi[1]
            t.set_up(energy_of_state_before_excitation=0, energy_of_state_after_excitation=1,
                     energy_zero_p = 0, energy_zero_s = 0,
                     orbital_to_excite_of=mo_p, orbital_to_excite_to=mo_s)
            assert self.m.symmetry_allowed_transition(transition=t) is combi[2]

