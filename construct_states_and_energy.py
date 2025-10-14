import itertools

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


def get_ground_state(norb,nel):
    if nel > 2*norb:
        raise Exception("too many electrons")
    state = []
    for i in range(norb):
        if nel >=2:
            state.append(2)
            nel -= 2
        elif nel == 0:
            state.append(0)
        else:
            state.append(1)
            nel -= 1
    return state



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


if __name__ == "__main__":
    x = construct_occupied_triplett_states(6,6)
    print(x)
