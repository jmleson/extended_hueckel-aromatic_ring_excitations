from TransitionIntegral import TransitionIntegral


class Transition:

    def __init__(self, n:int, s_after_transition:tuple[int,...], p_before_transition:tuple[int,...]):
        self.n = n
        if len(s_after_transition) != n or len(p_before_transition) != n:
            raise Exception("kjewrh")
        self.s_occupation_after_transition = s_after_transition
        self.p_before_transition = p_before_transition

        self.orbital_to_excite_of = None
        self.orbital_to_excite_to = None
        self.transition_integral = None

    def set_up(self, energy_of_state_before_excitation, energy_of_state_after_excitation,
               energy_zero_p, energy_zero_s,
               orbital_to_excite_of, orbital_to_excite_to):
        self.energy_of_state_before_excitation = energy_of_state_before_excitation
        self.energy_of_state_after_excitation = energy_of_state_after_excitation
        self.orbital_to_excite_to = orbital_to_excite_to
        self.orbital_to_excite_of = orbital_to_excite_of
        self.energy_zero_p = energy_zero_p
        self.energy_zeros_s = energy_zero_s
        self.transition_integral = TransitionIntegral(bra=self.orbital_to_excite_of, ket=self.orbital_to_excite_to)

    def get_difference_between_mean_orbital_energies(self):
        return self.energy_zeros_s - self.energy_zero_p

    def get_changed_orbital_index(self):
        if sum(self.s_occupation_after_transition) != 1:
            raise Exception("get_changed_orbital_index is only able to handle one excited electron")
        return [i for i, occ in enumerate(self.s_occupation_after_transition) if occ != 0][0]

    def get_s_before_transition(self):
        return tuple([0 for i in range(self.n)])

    def get_p_after_transition(self):
        return tuple([self.p_before_transition[i] - self.s_occupation_after_transition[i] for i in range(self.n)])

    def get_transitioning_energy(self):
        e1 = self.orbital_to_excite_to.eigenvalue - self.orbital_to_excite_of.eigenvalue# single occupation assumed
        e2 = self.energy_of_state_after_excitation - self.energy_of_state_before_excitation
        assert e1 == e2
        return e2

    def print(self):
        changed_orbital_index = self.get_changed_orbital_index()
        print(f"\t\033[1mtransition\033[0m from phi_{changed_orbital_index+1} into phi_s_{changed_orbital_index+1}:")
        print("\t\t-", self.p_before_transition, ",", self.get_s_before_transition(), "-->", self.get_p_after_transition(), ",", self.s_occupation_after_transition )
        if self.orbital_to_excite_of is not None:
                lc = 0
                for i in range(len(self.orbital_to_excite_of.salcs)):
                    lc += (self.orbital_to_excite_of.salcs[0]["factor"] * (self.orbital_to_excite_of.salcs[0]["salc"].equation ) )
                print(f"\t\t- phi_{changed_orbital_index+1}:\t linear combination =",lc,
                      f"\t\t, energy = {self.orbital_to_excite_of.eigenvalue}"
                      )
        if self.orbital_to_excite_to is not None:
                lc = 0
                for i in range(len(self.orbital_to_excite_to.salcs)):
                    lc += (self.orbital_to_excite_to.salcs[0]["factor"] * (self.orbital_to_excite_to.salcs[0]["salc"].equation))
                print(f"\t\t- xi_{changed_orbital_index + 1}:\t linear combination =", lc,
                      f"\t\t, energy = {self.orbital_to_excite_to.eigenvalue}"
                      )
        print("\t\t- Δ Energy:\t", self.get_transitioning_energy())
        print("\t\t- Energy of State before Excitation:\t", self.energy_of_state_before_excitation)
        print("\t\t- Energy of State after Excitation:\t", self.energy_of_state_after_excitation)
        print("\t\t- Energy between Average Energies of phi-/xi-orbitals:\t", self.get_difference_between_mean_orbital_energies())
        if self.transition_integral is not None:
            print("\t\t- dipole transition moment:\t", self.transition_integral.multiply_out())
        print()


import sympy as sp
def calculate_dispersion_energy(transitions_a:list[Transition], transitions_b:list[Transition]):
        r = sp.Symbol("r")
        T_zz = 2 / r ** (3)

        E_dispersion = 0
        for ia in transitions_a:
            for jb in transitions_b:
                zaehler = (ia.transition_integral.multiply_out() * T_zz * jb.transition_integral.multiply_out()) ** 2
                nenner = ia.get_transitioning_energy() + jb.get_transitioning_energy() # / (E_i^a-E_0 + E_j^b-E0)
                E_dispersion += zaehler/nenner
        return E_dispersion