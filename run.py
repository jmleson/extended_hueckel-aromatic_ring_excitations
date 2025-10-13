from MoleculeRepresentation import MoleculeRepresentation
from MoleculeState import MoleculeState



# m = MoleculeState(n=4)
# m.set_occupation(p_occupation=(2, 1, 1, 0))
# m.calculate_result_for_all_transitions_of_set_occupation()


# m.set_occupation(p_occupation=(2, 2, 0, 0))
# m.calculate_result_for_all_transitions_of_set_occupation()
#
# m.set_occupation(p_occupation=(2, 0, 2, 0))
# m.calculate_result_for_all_transitions_of_set_occupation()


m = MoleculeState(n=6)
m.set_occupation(p_occupation=(2,1,1,1,1,0))
# for i in m.s.get_energy_levels():
#     i.print()
m.calculate_result_for_all_transitions_of_set_occupation()




