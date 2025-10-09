from TransitionIntegral import TransitionIntegral


class Transition:

    def __init__(self, n:int, s_after_transition:tuple[int,...], p_before_transition:tuple[int,...]):
        self.n = n
        if len(s_after_transition) != n or len(p_before_transition) != n:
            raise Exception("kjewrh")
        self.s_occupation_after_transition = s_after_transition
        self.p_before_transition = p_before_transition

        self.transition_energy = None
        self.orbital_to_excite_of = None
        self.orbital_to_excite_to = None
        self.transition_integral = None

    def set_up(self, transition_energy, orbital_to_excite_of, orbital_to_excite_to):
        self.transition_energy = transition_energy
        self.orbital_to_excite_to = orbital_to_excite_to
        self.orbital_to_excite_of = orbital_to_excite_of
        self.transition_integral = TransitionIntegral(bra=self.orbital_to_excite_of, ket=self.orbital_to_excite_to)

    def get_changed_orbital_index(self):
        return [i for i in range(len(self.s_occupation_after_transition)) if i != 0][0]

    def get_s_before_transition(self):
        return tuple([0 for i in range(self.n)])

    def get_p_after_transition(self):
        return tuple([self.p_before_transition[i] - self.s_occupation_after_transition[i] for i in range(self.n)])


    def print(self):
        changed_orbital_index = self.get_changed_orbital_index()
        print(f"\033[1mtransition\033[0m from phi_{changed_orbital_index+1} into phi_s_{changed_orbital_index+1}:")
        print("\t- ", self.p_before_transition, ",", self.get_s_before_transition(), "-->", self.get_p_after_transition(), ",", self.s_occupation_after_transition )
        if self.orbital_to_excite_of is not None:
                lc = 0
                for i in range(len(self.orbital_to_excite_of.salcs)):
                    lc += (self.orbital_to_excite_of.salcs[0]["factor"] * (self.orbital_to_excite_of.salcs[0]["salc"].equation ) )
                print(f"\t- phi_{changed_orbital_index+1}:\t linear combination =",lc,
                      f"\t, energy = {self.orbital_to_excite_of.eigenvalue}"
                      )
        if self.orbital_to_excite_to is not None:
                lc = 0
                for i in range(len(self.orbital_to_excite_to.salcs)):
                    lc += (self.orbital_to_excite_to.salcs[0]["factor"] * (self.orbital_to_excite_to.salcs[0]["salc"].equation))
                print(f"\t- sigma_{changed_orbital_index + 1}:\t linear combination =", lc,
                      f"\t, energy = {self.orbital_to_excite_to.eigenvalue}"
                      )
        print("\t- Transition Energy:\t", self.transition_energy)
        if self.transition_integral is not None:
            print("\t- dipole transition moment:\t", self.transition_integral.multiply_out())
        print()

