import itertools

from MoleculeState import MoleculeState


# import sympy as sp

# from tst.solve_saekular_equation import get_six_orbitals_in_benzene_form


# def construct_states_and_energy(norb, nel):
#     states = construct_occupied_triplett_states(norb=norb, nel=nel)
#     p_orbitals = get_six_orbitals_in_benzene_form(alpha_symbol="alpha", beta_symbol="beta", info="p orbitals")#type: list[ molecule_orbital ]
#     p_orbital_energies = [i.eigenvalue for i in p_orbitals]
#
#     for state in states:
#         # print("state", state, end=" = ")
#         total_energy_of_state = sp.Integer(0)
#         for j in range(len(state)):
#             # if state[j] != 0:
#             #     print(f"+ {state[j]} * ({p_orbital_energies[j]})", end =" " )
#             total_energy_of_state += state[j] * (p_orbital_energies[j])
#         # print("=", total_energy_of_state)






# def calculate_result_for_all_transitions_print(triplet_states:list[dict], m:MoleculeState):
#     for triplett in triplet_states:
#         print(f"\033[1mState: {triplett}\033[0m")
#         m.set_occupation(p_occupation=triplett["occupation"])
#         allowed_transitions = m.calculate_result_for_all_transitions_of_set_occupation(print_active=True)
#


# if __name__ == "__main__":
#     x = construct_occupied_triplett_states(6,6)
#     print(x)
