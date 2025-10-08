import itertools
import sympy as sp

from IrreducibleRepresentation import IrreducibleRepresentation
from PointGroup import PointGroup
from SALC import norm_and_group_SALCs
from SymmetryOperation import SymmetryOperation
from tst.solve_saekular_equation import calculate

###### SET UP TESTING CASE ###############
p = PointGroup(n=4)


##############################
reducible_representation = p.get_reducible_representation_for_ring_p_orbitals()
print(reducible_representation)

# irreducible_representation = p.decomposing_into_irreducible_representations(reducible_representation)
# print(irreducible_representation)


# p.get_energy_levels()








