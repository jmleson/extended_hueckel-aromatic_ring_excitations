from IrreducibleRepresentation import IrreducibleRepresentation
from PointGroups.PointGroup import PointGroup
from SymmetryOperation import SymmetryOperation


class C2v(PointGroup):
    def __init__(self, n:int=3):
        super().__init__(order=16, n=4)
        self.set_up_symmetry_operations()
        self.set_up_irreducible_representations()
        self.total_symmetric_representation = "A1"

    def set_up_symmetry_operations(self):
        if self.n == 4:
            self.operations = []
            o = SymmetryOperation(n=self.n, name="E", transform_p=lambda i: i, amount=1)
            self.operations.append(o)
            o = SymmetryOperation(n=self.n, name="C2", transform_p=lambda i: [-4, -3, -2, -1][i - 1], amount=1)
            self.operations.append(o)
            o = SymmetryOperation(n=self.n, name="σv(xz)", transform_p=lambda i: [-1, -2, -3, -4][i - 1], amount=1)
            self.operations.append(o)
            o = SymmetryOperation(n=self.n, name="σv(yz)", transform_p=lambda i: [4, 3, 2, 1][i - 1], amount=1)
            self.operations.append(o)

            self.dipole_operator_symmetry = "B2"# 1  -1   -1  1

    def set_up_irreducible_representations(self):
        if self.n == 4:
            self.irreducible_representations = []
            symmetry_names = [op.name for op in self.operations]
            i = IrreducibleRepresentation({j: 1 for j in symmetry_names}, name="A1")
            self.irreducible_representations.append(i)
            A2_chars = [1, 1, -1, -1]
            i = IrreducibleRepresentation(dict(zip(symmetry_names, A2_chars)), name="A2")
            self.irreducible_representations.append(i)
            B1_chars = [1, -1, 1, -1]
            i = IrreducibleRepresentation(dict(zip(symmetry_names, B1_chars)), name="B1")
            self.irreducible_representations.append(i)
            B2_chars = [1, -1, -1, 1]
            i = IrreducibleRepresentation(dict(zip(symmetry_names, B2_chars)), name="B2")
            self.irreducible_representations.append(i)
