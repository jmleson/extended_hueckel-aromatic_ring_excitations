import itertools

from get_ground_state import get_ground_state


def construct_occupied_triplett_states(norb, nel):
    """
    Generate all possible occupation states for `norb` orbitals and `nel` electrons.
    Each orbital can have 0, 1, or 2 electrons.
    Only states with at least 2 unpaired electrons (two 1s) are included.
    """
    all_combinations = itertools.product([0, 1, 2], repeat=norb-1)

    valid_states = [get_ground_state(norb=norb, nel=nel)]
    for state in all_combinations:
        total_electrons = sum(state)
        unpaired_electrons = state.count(1)
        excited = sum(state[3:])

        if total_electrons == nel and unpaired_electrons >= 2 and (
                    (excited <= 1 or (excited == 2 and state[3] == 1 and state[4] == 1 and state[2] == 1 and state[1] == 1))
                    and state[0] == 2
        ):
            # Optional: Auffüllen auf feste Länge, falls benötigt
            # padded_state = state + (0,)*(3-len(state))
            valid_states.append( state+ (0,) )


    return [{"occupation": state,
             "unpaired electrons": state.count(1),
             "multiplicity": state.count(1)*(1/2)*2 + 1 # assumes alpha-spin for all unpaired electrons
            } for state in valid_states]