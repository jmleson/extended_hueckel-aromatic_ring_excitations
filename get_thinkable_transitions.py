import copy


def get_thinkable_transitions(state):
    """
    Generate all possible single-electron transitions from a given state.
    Each transition is a tuple (from_state, to_state), where:
    - from_state: state after removing one electron from an orbital
    - to_state: state after adding one electron to another orbital
    """
    transitions = []
    for i in range(1, len(state)):# keine Anregungen aus unterstem besetzten Orbital
        if state[i] > 0:
            from_state = list(copy.deepcopy(state))
            from_state[i] -= 1
            for j in range(len(state)):
                    to_state = [0 for i in range(len(state))]
                    to_state[j] += 1
                    transition = (state, "->", from_state, to_state)
                    transitions.append(transition)
    return transitions




if __name__ == '__main__':
    state = (2, 2, 1, 1, 0, 0)

    thinkable_transitions = get_thinkable_transitions(state)
    for t in thinkable_transitions:
        print(t)

    assert len(thinkable_transitions) == 3*6
