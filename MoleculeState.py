
from Transition import Transition
from molecule_orbital import molecule_orbital
import sympy as sp


class MoleculeState:
    def __init__(self, bonding_p: list[molecule_orbital], antibonding_s: list[molecule_orbital]):
        self.bonding_p = bonding_p
        self.n = len(self.bonding_p)

        alpha, alpha_s, beta, beta_s = sp.symbols("alpha alpha_s beta beta_s")
        self.antibonding_s = []
        for s in range(len(antibonding_s)):
            antibonding_s[s].eigenvalue = antibonding_s[s].eigenvalue.subs(alpha, alpha_s)
            antibonding_s[s].eigenvector = antibonding_s[s].eigenvector.subs(beta, beta_s)
            self.antibonding_s.append(antibonding_s[s])

        self.p_occupation = None

    def set_occupation(self, p_occupation):
        self.p_occupation = p_occupation

    def calculate_energy_before_transition(self):
        if self.p_occupation is None:
            raise Exception("perform set_occupation() before")
        return self.calculate_energy_for_state_and_occupation(state=self.bonding_p, occupation=self.p_occupation)

    def calculate_energy_for_state_and_occupation(self, state:list[molecule_orbital], occupation:tuple[int,...]):
        if len(state) != self.n or len(occupation) != self.n:
            raise Exception("state and occupation must have same length")
        e = 0
        for i in range(self.n):
            e += state[i].eigenvalue * occupation[i]
        return e


    def get_p_after_transition(self, s_after_transition:tuple[int]):
        return tuple([self.p_occupation[i] - s_after_transition[i] for i in range(len(self.p_occupation))])


    def get_thinkable_transitions(self) -> list[Transition]:
        if self.p_occupation is None:
            raise Exception("perform set_occupation() before")
        if len(self.p_occupation) != self.n:
            raise Exception(f"unfitting occupation sequence for {self.n} orbitals")
        e_state_before = self.calculate_energy_before_transition()
        # 1st orbital (lowest energy) unchanged
        s_occupations = []
        for i in range(1, len(self.p_occupation)):
            if self.p_occupation[i] > 0:
                s_after_transition = tuple([ 1 if p == i else 0 for p in range(len(self.p_occupation)) ])
                transition = Transition(n=self.n, s_after_transition=s_after_transition,
                                        p_before_transition=self.p_occupation)

                E_s_state = self.calculate_energy_for_state_and_occupation(state=self.antibonding_s,
                                                                   occupation=s_after_transition)
                energy_single_p = self.calculate_energy_for_state_and_occupation(state=self.bonding_p,
                                                                   occupation=s_after_transition)
                E_state_after = E_s_state + e_state_before - energy_single_p
                changed_orbital_index = transition.get_changed_orbital_index()
                transition.set_up(orbital_to_excite_of = self.bonding_p[changed_orbital_index],
                                  orbital_to_excite_to = self.antibonding_s[changed_orbital_index],
                                  energy_of_state_before_excitation = e_state_before,
                                  energy_of_state_after_excitation = E_state_after)

                assert E_s_state - energy_single_p == E_state_after - e_state_before
                assert transition.get_transitioning_energy() == E_s_state - energy_single_p

                s_occupations.append(transition)
        return s_occupations


    def symmetry_allowed_transition(self):
        pass#TODO

    def calculate_result_for_all_transitions_of_set_occupation(self):
        transitions = self.get_thinkable_transitions()

        for transition in transitions:
            transition.print()