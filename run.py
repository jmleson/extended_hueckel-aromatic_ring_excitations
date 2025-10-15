import sympy as sp

from MoleculeState import MoleculeState
from Transition import calculate_dispersion_energy
from construct_states_and_energy import construct_occupied_triplett_states

# m = MoleculeState(n=4)
# m.set_occupation(p_occupation=(2, 1, 1, 0))
# m.calculate_result_for_all_transitions_of_set_occupation()


# m.set_occupation(p_occupation=(2, 2, 0, 0))
# m.calculate_result_for_all_transitions_of_set_occupation()
#
# m.set_occupation(p_occupation=(2, 0, 2, 0))
# m.calculate_result_for_all_transitions_of_set_occupation()


m = MoleculeState(n=6)


triplet_states = construct_occupied_triplett_states(nel=6, norb=6)
for triplett in triplet_states:
    print(f"\033[1mState: {triplett}\033[0m")
    m.set_occupation(p_occupation=triplett["occupation"])
    m.calculate_result_for_all_transitions_of_set_occupation(print_active=True)


state_no = 0
for a in triplet_states:
    for b in triplet_states:
        if a["multiplicity"] + b["multiplicity"] == 5+1:
            state_no += 1
            print(f"{state_no}) A =", a["occupation"], "<-> B=", b["occupation"] )
            m.set_occupation(p_occupation=a["occupation"])
            triplett_states_A = m.calculate_result_for_all_transitions_of_set_occupation(print_active=False)
            m.set_occupation(p_occupation=b["occupation"])
            triplett_states_B = m.calculate_result_for_all_transitions_of_set_occupation(print_active=False)
            e_disp = calculate_dispersion_energy(triplett_states_A, triplett_states_B)
            print("\tE_dispersion =", e_disp)
        else:
            # skip triplet quintet combinations
            pass
            # print("Multiplicity",a["multiplicity"] + b["multiplicity"], ", A =", a["occupation"], "<-> B=", b["occupation"] )












