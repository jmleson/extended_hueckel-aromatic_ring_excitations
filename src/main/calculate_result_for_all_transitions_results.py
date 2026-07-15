from MoleculeState import MoleculeState
from Transition import calculate_dispersion_energy


def calculate_result_for_all_transitions_results(m:MoleculeState, triplet_states:list[dict]):
    state_no = 0
    for a in triplet_states:
        for b in triplet_states:
            if a["multiplicity"] + b["multiplicity"] == 5 + 1:
                state_no += 1
                print(f"{state_no}) A =", a["occupation"], "<-> B =", b["occupation"])
                m.set_occupation(p_occupation=a["occupation"])
                triplett_states_A_with_warning = m.calculate_result_for_all_transitions_of_set_occupation(print_active=False)
                triplett_states_A = [t[0] for t in triplett_states_A_with_warning]
                m.set_occupation(p_occupation=b["occupation"])
                triplett_states_B_with_warning = m.calculate_result_for_all_transitions_of_set_occupation(print_active=False)
                triplett_states_B = [t[0] for t in triplett_states_B_with_warning]
                e_disp = calculate_dispersion_energy(triplett_states_A, triplett_states_B)
                print("\tE_dispersion =", e_disp)
            else:
                # skip triplet quintet combinations
                pass
                # print("Multiplicity",a["multiplicity"] + b["multiplicity"], ", A =", a["occupation"], "<-> B=", b["occupation"] )