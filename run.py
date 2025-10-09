
from MoleculeRepresentation import MoleculeRepresentation
from MoleculeState import MoleculeState
from TransitionIntegral import TransitionIntegral
from molecule_orbital import molecule_orbital
import sympy as sp


p = MoleculeRepresentation(n=4)


# print("P ORBITALS")
p_mo_orbitals = p.get_energy_levels()

# for x in p_mo_orbitals:
#     print("\t", x.symmetry, ":\t", x.eigenvalue)


# print("\nS ORBITALS")
p.set_to_s_orbitals()
s_mo_orbitals = p.get_energy_levels()

# for s in s_mo_orbitals:
#     print("\t", s.symmetry, ":\t", s.eigenvalue)





# for px in p_mo_orbitals:
#     px.print()


m = MoleculeState(bonding_p=p_mo_orbitals, antibonding_s=s_mo_orbitals)
m.set_occupation(p_occupation=(2, 0, 1, 0)) # phi_3
# print("E before transition", m.calculate_energy_before_transition())


m.calculate_result_for_all_transitions_of_set_occupation()






