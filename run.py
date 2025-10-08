
from MoleculeRepresentation import MoleculeRepresentation
from molecule_orbital import molecule_orbital
import sympy as sp


p = MoleculeRepresentation(n=4)

print("P ORBITALS")
p_mo_orbitals = p.get_energy_levels()

for x in p_mo_orbitals:
    print("\t", x.symmetry, ":\t", x.eigenvalue)


print("\nS ORBITALS")
p.set_to_s_orbitals()
s_mo_orbitals = p.get_energy_levels()

for s in s_mo_orbitals:
    print("\t", s.symmetry, ":\t", s.eigenvalue)




class MoleculeState:
    def __init__(self, bonding_p: list[molecule_orbital], antibonding_s: list[molecule_orbital]):
        self.bonding_p = bonding_p
        self.n = len(self.bonding_p)

        alpha, alpha_s, beta, beta_s = sp.symbols("alpha alpha_s beta beta_s")
        self.antibonding_s = []
        for s in range(len(antibonding_s)):
            antibonding_s[s].eigenvalue = antibonding_s[s].eigenvalue.subs(alpha, alpha_s)
            antibonding_s[s].eigenvector = antibonding_s[s].eigenvector.subs(alpha, alpha_s)
            self.antibonding_s.append(antibonding_s[s])

        self.p_occupation = None

    def set_occupation(self, p_occupation):
        self.p_occupation = p_occupation

    def calculate_energy_before_transition(self):
        if self.p_occupation is None:
            raise Exception("perform set_occupation() before")
        return self.calculate_energy_for_state_and_occupation(state=self.bonding_p, occupation=self.p_occupation)

    def calculate_energy_for_state_and_occupation(self, state:list[molecule_orbital], occupation:tuple[int]):
        if len(state) != self.n or len(occupation) != self.n:
            raise Exception("state and occupation must have same length")
        e = 0
        for i in range(self.n):
            e += state[i].eigenvalue * occupation[i]
        return e


    def get_p_after_transition(self, s_after_transition:tuple[int]):
        return tuple([self.p_occupation[i] - s_after_transition[i] for i in range(len(self.p_occupation))])


    def get_thinkable_transitions(self):
        if self.p_occupation is None:
            raise Exception("perform set_occupation() before")
        if len(self.p_occupation) != self.n:
            raise Exception(f"unfitting occupation sequence for {self.n} orbitals")
        e_before = self.calculate_energy_before_transition()
        # 1st orbital (lowest energy) unchanged
        s_occupations = []
        for i in range(1, len(self.p_occupation)):
            if self.p_occupation[i] > 0:
                s_after_transition = tuple([ 1 if p == i else 0 for p in range(len(self.p_occupation)) ])
                E = self.calculate_energy_for_state_and_occupation(state=self.antibonding_s, occupation=s_after_transition)
                s_occupations.append({"s occupation after transition": s_after_transition,
                                      "Delta E by transition": E - e_before
                                      })
        return s_occupations


    def symmetry_allowed_transition(self):
        pass







m = MoleculeState(bonding_p=p_mo_orbitals, antibonding_s=s_mo_orbitals)
m.set_occupation(p_occupation=(2, 1, 1, 0))
print(m.calculate_energy_before_transition())

print( m.get_thinkable_transitions() )






