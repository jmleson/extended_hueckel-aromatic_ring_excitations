


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

