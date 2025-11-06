from IrreducibleRepresentation import IrreducibleRepresentation
from PointGroups.PointGroup import PointGroup
from SymmetryOperation import SymmetryOperation


class C2v(PointGroup):
    def __init__(self, n:int=3):
        super().__init__(order=4, n=n)

        self.set_up_symmetry_operations()
        self.set_up_irreducible_representations()
        self.total_symmetric_representation = "A1"

    def set_up_symmetry_operations(self):
        self.operations = []
        o = SymmetryOperation(n=self.n, name="E", transform_p=lambda i: i, amount=1)
        self.operations.append(o)
        if self.n == 4:
            o = SymmetryOperation(n=self.n, name="C2", transform_p=lambda i: [-4, -3, -2, -1][i - 1], amount=1)
            self.operations.append(o)
            o = SymmetryOperation(n=self.n, name="σv(xz)", transform_p=lambda i: [-1, -2, -3, -4][i - 1], amount=1)
            self.operations.append(o)
            o = SymmetryOperation(n=self.n, name="σv(yz)", transform_p=lambda i: [4, 3, 2, 1][i - 1], amount=1)
            self.operations.append(o)

            self.dipole_operator_symmetry = "B2"# 1  -1   -1  1
        elif self.n == 7:#Chlorobenzene
            # last item == Cl
            o = SymmetryOperation(n=self.n, name="C2", transform_p=lambda i: [-1, -6, -5, -4, -3, -2, -7][i - 1], amount=1)
            self.operations.append(o)
            o = SymmetryOperation(n=self.n, name="σv(xz)", transform_p=lambda i: [-1, -2, -3, -4, -5, -6, -7][i - 1], amount=1)
            self.operations.append(o)
            o = SymmetryOperation(n=self.n, name="σv(yz)", transform_p=lambda i: [1, 6, 5, 4, 3, 2, 7][i - 1], amount=1)
            self.operations.append(o)

            self.dipole_operator_symmetry = "B2"  #1 -1 -1 1
        elif self.n == 6: # Pyridin
            # last item == Cl
            o = SymmetryOperation(n=self.n, name="C2", transform_p=lambda i: [-1, -6, -5, -4, -3, -2][i - 1], amount=1)
            self.operations.append(o)
            o = SymmetryOperation(n=self.n, name="σv(xz)", transform_p=lambda i: [-1, -2, -3, -4, -5, -6][i - 1], amount=1)
            self.operations.append(o)
            o = SymmetryOperation(n=self.n, name="σv(yz)", transform_p=lambda i: [1, 6, 5, 4, 3, 2][i - 1], amount=1)
            self.operations.append(o)

            self.dipole_operator_symmetry = "B2"  # 1 -1 -1 1
        else:
            raise Exception("point group case not yet implemented")

    def set_up_irreducible_representations(self):
        # if self.n == 4:
            self.irreducible_representations = []
            symmetry_names = list(dict.fromkeys(op.name for op in self.operations))
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

    def multiply(self, irred_1: str, irred_2: str) -> list[str]:
        possible_irreducible_representations = [i.name for i in self.irreducible_representations]
        if irred_1 not in possible_irreducible_representations or irred_2 not in possible_irreducible_representations:
            raise Exception(f"unknown irreducible representation {irred_1} or {irred_2}")

        if irred_1 == irred_2:
            return ["A1"]
        if irred_1 == "A2":
            return [self.replace_number(irred_2)]
        if irred_2 == "A2":
            return [self.replace_number(irred_1)]

        if irred_1 == "B1":
            return [self.replace_symbol(irred_2)]
        if irred_2 == "B1":
            return [self.replace_symbol(irred_1)]

        if irred_1 == "B2":
            return [self.replace_symbol(self.replace_number(irred_2))]
        if irred_2 == "B2":
            return [self.replace_symbol(self.replace_number(irred_1))]

if __name__ == "__main__":
    C2v = C2v(n=7)

    C2v.print_multiplication_table()