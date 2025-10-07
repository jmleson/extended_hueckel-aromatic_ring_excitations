from IrreducibleRepresentation import IrreducibleRepresentation
from PointGroup import PointGroup
from SymmetryOperation import SymmetryOperation

###### SET UP TESTING CASE ###############
p = PointGroup(n=6)
# revert Point group to C2v for 1,3-Butadien:
p.n = 4
# set_up_symmetry_operations:
p.operations = []
o = SymmetryOperation(n=p.n, name="E", transform_p=lambda i: i, amount=1)
p.operations.append(o)
o = SymmetryOperation(n=p.n, name="C2", transform_p=lambda i: [-4, -3, -2, -1][i - 1], amount=1)
p.operations.append(o)
o = SymmetryOperation(n=p.n, name="σv(xz)", transform_p=lambda i: [-1, -2, -3, -4][i - 1], amount=1)
p.operations.append(o)
o = SymmetryOperation(n=p.n, name="σv(yz)", transform_p=lambda i: [4, 3, 2, 1][i - 1], amount=1)
p.operations.append(o)

#set_up_irreducible_representations:
p.irreducible_representations = []
symmetry_names = [op.name for op in p.operations]
i = IrreducibleRepresentation({j: 1 for j in symmetry_names}, name="A1")
p.irreducible_representations.append(i)
A2_chars = [1, 1, -1, -1]
i = IrreducibleRepresentation(dict(zip(symmetry_names, A2_chars)), name="A2")
p.irreducible_representations.append(i)
B1_chars = [1, -1, 1, -1]
i = IrreducibleRepresentation(dict(zip(symmetry_names, B1_chars)), name="B1")
p.irreducible_representations.append(i)
B2_chars = [1, -1, -1, 1]
i = IrreducibleRepresentation(dict(zip(symmetry_names, B2_chars)), name="B2")
p.irreducible_representations.append(i)


##############################
# reducible_representation = p.get_reducible_representation_for_ring_p_orbitals()
# print(reducible_representation)
#
# irreducible_representation = p.decomposing_into_irreducible_representations(reducible_representation)
# print(irreducible_representation)

SALCs = p.get_all_SALCs()



