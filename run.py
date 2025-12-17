
from MoleculeState import MoleculeState


# def cyclopentadiene():
#     m = MoleculeState(n=4, bound_cl_to_c_positions=[], n_instead_of_c=[])
#     m.set_occupation(p_occupation=(2, 1, 1, 0))
#     m.calculate_result_for_all_transitions_of_set_occupation()
#
#     m.set_occupation(p_occupation=(2, 2, 0, 0))
#     m.calculate_result_for_all_transitions_of_set_occupation()
#
#     m.set_occupation(p_occupation=(2, 0, 2, 0))
#     m.calculate_result_for_all_transitions_of_set_occupation()



def benzene(print_active:bool=True):
    m = MoleculeState(n=6, bound_cl_to_c_positions=[],n_instead_of_c=[])
    if not print_active:
        m.latex_datei_erstellen("Benzene")
    else:
        # m.calculate_result_for_all_transitions(print_active=print_active)
        # m.calculate_result_for_all_transitions_results(print_active=print_active)
        m.compare_ground_state_assumption(print_active=print_active)



def pyridine(print_active:bool=True):
    m = MoleculeState(n=6, bound_cl_to_c_positions=[], n_instead_of_c=[1] )
    if not print_active:
        m.latex_datei_erstellen("Pyridine")
    else:
        m.calculate_result_for_all_transitions(print_active=print_active)
        m.calculate_result_for_all_transitions_results(print_active=print_active)
        m.compare_ground_state_assumption(print_active=print_active)

def pyrazine(print_active:bool=True):
    m = MoleculeState(n=6, bound_cl_to_c_positions=[], n_instead_of_c=[1, 4] )
    if not print_active:
        m.latex_datei_erstellen("Pyrazine")
    else:
        m.calculate_result_for_all_transitions(print_active=print_active)
        m.calculate_result_for_all_transitions_results(print_active=print_active)
        m.compare_ground_state_assumption(print_active=print_active)


def chlorobenzene(print_active:bool=True):
    m = MoleculeState(n=6, bound_cl_to_c_positions=[1], n_instead_of_c=[])
    if not print_active:
        m.latex_datei_erstellen("Chlorobenzene")
    else:
        m.calculate_result_for_all_transitions(print_active=print_active)
        m.calculate_result_for_all_transitions_results(print_active=print_active)
        m.compare_ground_state_assumption(print_active=print_active)

def dichlorobenzene(print_active:bool=True):
    m = MoleculeState(n = 6, bound_cl_to_c_positions=[1,4], n_instead_of_c=[])
    if not print_active:
        m.latex_datei_erstellen("Dichlorobenzene")
    else:
        m.calculate_result_for_all_transitions(print_active=print_active)
        m.calculate_result_for_all_transitions_results(print_active=print_active)
        m.compare_ground_state_assumption(print_active=print_active)

def hexachlorobenzene(print_active:bool=True):
    m = MoleculeState(n = 6, bound_cl_to_c_positions=[1,2,3,4,5,6], n_instead_of_c=[])
    if not print_active:
        m.latex_datei_erstellen("Hexachlorobenzene")
    else:
        m.calculate_result_for_all_transitions(print_active=print_active)
        m.calculate_result_for_all_transitions_results(print_active=print_active)
        m.compare_ground_state_assumption(print_active=print_active)

def hexafluorobenzene(print_active:bool=True):
    m = MoleculeState(n = 6, bound_cl_to_c_positions=[1,2,3,4,5,6], n_instead_of_c=[], heterosymbol = "F")
    if not print_active:
        m.latex_datei_erstellen("Hexafluorobenzene")
    else:
        m.calculate_result_for_all_transitions(print_active=print_active)
        m.calculate_result_for_all_transitions_results(print_active=print_active)
        m.compare_ground_state_assumption(print_active=print_active)



if __name__ == "__main__":
    print_active = False

    # cyclopentadiene()
    benzene(print_active)
    chlorobenzene(print_active)
    pyridine(print_active)
    pyrazine(print_active)

    dichlorobenzene(print_active)
    hexachlorobenzene(print_active)
    hexafluorobenzene(print_active)
