
from MoleculeState import MoleculeState
from construct_states_and_energy import construct_occupied_triplett_states, get_ground_state, \
    calculate_result_for_all_transitions_print, calculate_result_for_all_transitions_results, \
    compare_ground_state_assumption


def cyclopentadiene():
    m = MoleculeState(n=4, bound_cl_to_c_positions=[], n_instead_of_c=[])
    m.set_occupation(p_occupation=(2, 1, 1, 0))
    m.calculate_result_for_all_transitions_of_set_occupation()

    m.set_occupation(p_occupation=(2, 2, 0, 0))
    m.calculate_result_for_all_transitions_of_set_occupation()

    m.set_occupation(p_occupation=(2, 0, 2, 0))
    m.calculate_result_for_all_transitions_of_set_occupation()



def benzene():
    m = MoleculeState(n=6, bound_cl_to_c_positions=[],n_instead_of_c=[])
    triplet_states = construct_occupied_triplett_states(nel=6, norb=6)
    calculate_result_for_all_transitions_print(triplet_states=triplet_states)
    calculate_result_for_all_transitions_results(m=m, triplet_states=triplet_states)

    ground_state_occupation = get_ground_state(nel=6, norb=6)
    compare_ground_state_assumption(m=m, ground_state_occupation=ground_state_occupation)



def pyridine():
    m = MoleculeState(n=6, bound_cl_to_c_positions=[], n_instead_of_c=[1] )
    triplet_states = construct_occupied_triplett_states(nel=6, norb=6)
    calculate_result_for_all_transitions_print(triplet_states=triplet_states, m=m)
    calculate_result_for_all_transitions_results(m=m, triplet_states=triplet_states)

    ground_state_occupation = get_ground_state(nel=6, norb=6)
    compare_ground_state_assumption(m=m, ground_state_occupation=ground_state_occupation)

def pyrazine():
    m = MoleculeState(n=6, bound_cl_to_c_positions=[], n_instead_of_c=[1, 3] )
    triplet_states = construct_occupied_triplett_states(nel=6, norb=6)
    calculate_result_for_all_transitions_print(triplet_states=triplet_states, m=m)
    calculate_result_for_all_transitions_results(m=m, triplet_states=triplet_states)

    ground_state_occupation = get_ground_state(nel=6, norb=6)
    compare_ground_state_assumption(m=m, ground_state_occupation=ground_state_occupation)


def chlorobenzene():
    m = MoleculeState(n=6, bound_cl_to_c_positions=[1], n_instead_of_c=[])
    triplet_states = construct_occupied_triplett_states(nel=6, norb=6)
    # adjust to full Cl orbital:
    for triplett in triplet_states:
        triplett["occupation"] += (2,)

    calculate_result_for_all_transitions_print(triplet_states=triplet_states, m=m)
    calculate_result_for_all_transitions_results(m=m, triplet_states=triplet_states)
    ground_state_occupation = get_ground_state(nel=6, norb=6)+[2]
    compare_ground_state_assumption(m=m, ground_state_occupation=ground_state_occupation)



if __name__ == "__main__":
    # cyclopentadiene()
    # benzene()
    # chlorobenzene()
    # pyridine()
    pyrazine()