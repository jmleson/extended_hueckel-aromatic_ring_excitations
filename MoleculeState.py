from MoleculeRepresentation import MoleculeRepresentation
from Transition import Transition
from molecule_orbital import molecule_orbital
import sympy as sp


class MoleculeState:
    def __init__(self, n:int):
        self.n = n
        self.p = MoleculeRepresentation(n=n)

        self.s = MoleculeRepresentation(n=n)
        self.set_up()

        self.p_occupation = None

    def set_up(self):
        # print("P ORBITALS")
        p_mo_orbitals = self.p.get_energy_levels()
        # for x in p_mo_orbitals:
        #     print("\t", x.symmetry, ":\t", x.eigenvalue)
        self.bonding_p = p_mo_orbitals

        # print("\nS ORBITALS")
        self.s.set_to_s_orbitals()
        s_mo_orbitals = self.s.get_energy_levels()
        # for s in s_mo_orbitals:
        #     print("\t", s.symmetry, ":\t", s.eigenvalue)
        self.antibonding_s = []
        alpha, alpha_s, beta, beta_s = sp.symbols("alpha alpha_s beta beta_s")
        for s in range(len(s_mo_orbitals)):
            s_mo_orbitals[s].eigenvalue = s_mo_orbitals[s].eigenvalue.subs(alpha, alpha_s).subs(beta, beta_s)
            s_mo_orbitals[s].eigenvector = s_mo_orbitals[s].eigenvector.subs(alpha, alpha_s).subs(beta, beta_s)
            self.antibonding_s.append(s_mo_orbitals[s])

    def set_occupation(self, p_occupation):
        if len(p_occupation) != self.n:
            raise Exception("wrong number of orbitals")
        for n in p_occupation:
            if n not in [0,1,2]:
                raise Exception("wrong occupation")
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
            state[i].occupation = occupation[i]
            # state[i].print()
            e += state[i].getEnergy()
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
                E_state_after = E_s_state + (e_state_before - energy_single_p)
                changed_orbital_index = transition.get_changed_orbital_index()
                transition.set_up(orbital_to_excite_of = self.bonding_p[changed_orbital_index],
                                  orbital_to_excite_to = self.antibonding_s[changed_orbital_index],
                                  energy_zero_p = self.calculate_mean_energy_for_state(state=self.bonding_p),
                                  energy_zero_s = self.calculate_mean_energy_for_state(state=self.antibonding_s),
                                  energy_of_state_before_excitation = e_state_before,
                                  energy_of_state_after_excitation = E_state_after)

                assert E_s_state - energy_single_p == E_state_after - e_state_before
                assert transition.get_transitioning_energy() == E_state_after - e_state_before

                s_occupations.append(transition)
        return s_occupations

    def calculate_mean_energy_for_state(self, state):
        # calculate "zero" energy of orbital set, as defined in state (p_bonding / s_antibonding)
        energy_of_single_total_occupation = self.calculate_energy_for_state_and_occupation(state=state,
                                                    occupation=[1 for i in range(len(state))])
        return energy_of_single_total_occupation / len(state)

    def symmetry_allowed_transition(self, transition:Transition):

        sym_p = transition.orbital_to_excite_of.symmetry
        sym_s = transition.orbital_to_excite_to.symmetry
        dipole_operator = self.p.pointgroup.dipole_operator_symmetry
        # calculate symmetry of total integral: sym_p  x  dipole_operator  x  sym_s
        part_I = self.p.pointgroup.multiply(dipole_operator, sym_s)# new pseudo p orbital
        part_II = [self.p.pointgroup.multiply(sym_p, i) for i in part_I]
        # print(dipole_operator, "x", sym_s, "x", sym_p,  "=", part_II)
        return self.p.pointgroup.total_symmetric_representation in [x for row in part_II for x in row]

    def calculate_result_for_all_transitions_of_set_occupation(self, print_active:bool=True):
        transitions = self.get_thinkable_transitions()

        allowed_transitions = []
        for transition in transitions:
            if self.symmetry_allowed_transition(transition):
                if print_active:
                    transition.print()
                allowed_transitions.append(transition)
            else:
                if print_active:
                    print("\tFORBIDDEN\n")
            #     transition.print()
        return allowed_transitions


