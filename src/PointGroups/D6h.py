from src.main.IrreducibleRepresentation import IrreducibleRepresentation
from src.PointGroups.PointGroup import PointGroup
from src.main.SymmetryOperation import SymmetryOperation


class D6h(PointGroup):
    def __init__(self, n: int = 6):
        super().__init__(order=24, n=n)
        self.set_up_symmetry_operations()
        self.set_up_irreducible_representations()
        self.total_symmetric_representation = "A1g"

    def set_up_symmetry_operations(self):
        self.operations = []
        o = SymmetryOperation(n=self.n, name="E", transform_p=lambda i: i, amount=1)
        self.operations.append(o)
        if self.n == 12:
            o = SymmetryOperation(n=self.n, name="C6", transform_p=lambda i: [6, 1, 2, 3, 4, 5, 12, 7, 8, 9, 10, 11][i - 1], amount=1)
            self.operations.append(o)
            o = SymmetryOperation(n=self.n, name="C6", transform_p=lambda i: [2, 3, 4, 5, 6, 1, 8, 9, 10, 11, 12, 7][i - 1], amount=1)
            self.operations.append(o)

            o = SymmetryOperation(n=self.n, name="C3", transform_p=lambda i: [5, 6, 1, 2, 3, 4, 11, 12, 7, 8, 9, 10][i - 1], amount=1)
            self.operations.append(o)
            o = SymmetryOperation(n=self.n, name="C3", transform_p=lambda i: [3, 4, 5, 6, 1, 2, 9, 10, 11, 12, 7, 8][i - 1], amount=1)
            self.operations.append(o)

            o = SymmetryOperation(n=self.n, name="C2", transform_p=lambda i: [4, 5, 6, 1, 2, 3, 10, 11, 12, 7, 8, 9][i - 1], amount=1)
            self.operations.append(o)

            o = SymmetryOperation(n=self.n, name="C\'2", transform_p=lambda i: [-1, -6, -5, -4, -3, -2, -7, -12, -11, -10, -9, -8][i - 1],
                                  amount=1)
            self.operations.append(o)
            o = SymmetryOperation(n=self.n, name="C\'2", transform_p=lambda i: [-3, -2, -1, -6, -5, -4, -9, -8, -7, -12, -11, -10][i - 1],
                                  amount=1)
            self.operations.append(o)
            o = SymmetryOperation(n=self.n, name="C\'2", transform_p=lambda i: [-5, -4, -3, -2, -1, -6, -11, -10, -9, -8, -7, -12][i - 1],
                                  amount=1)
            self.operations.append(o)

            o = SymmetryOperation(n=self.n, name="C\'\'2", transform_p=lambda i: [-2, -1, -6, -5, -4, -3, -8, -7, -12, -11, -10, -9][i - 1],
                                  amount=1)
            self.operations.append(o)
            o = SymmetryOperation(n=self.n, name="C\'\'2", transform_p=lambda i: [-4, -3, -2, -1, -6, -5, -10, -9, -8, -7, -12, -11][i - 1],
                                  amount=1)
            self.operations.append(o)
            o = SymmetryOperation(n=self.n, name="C\'\'2", transform_p=lambda i: [-6, -5, -4, -3, -2, -1, -12, -11, -10, -9, -8, -7][i - 1],
                                  amount=1)
            self.operations.append(o)

            o = SymmetryOperation(n=self.n, name="i", transform_p=lambda i: [-4, -5, -6, -1, -2, -3, -10, -11, -12, -7, -8, -9][i - 1], amount=1)
            self.operations.append(o)

            o = SymmetryOperation(n=self.n, name="S3", transform_p=lambda i: [-5, -6, -1, -2, -3, -4, -11, -12, -7, -8, -9, -10][i - 1], amount=1)
            self.operations.append(o)
            o = SymmetryOperation(n=self.n, name="S3", transform_p=lambda i: [-3, -4, -5, -6, -1, -2, -9, -10, -11, -12, -7, -8][i - 1], amount=1)
            self.operations.append(o)

            o = SymmetryOperation(n=self.n, name="S6", transform_p=lambda i: [-6, -1, -2, -3, -4, -5, -12, -7, -8, -9, -10, -11][i - 1], amount=1)
            self.operations.append(o)
            o = SymmetryOperation(n=self.n, name="S6", transform_p=lambda i: [-2, -3, -4, -5, -6, -1, -8, -9, -10, -11, -12, -7][i - 1], amount=1)
            self.operations.append(o)

            o = SymmetryOperation(n=self.n, name="σh", transform_p=lambda i: [-1, -2, -3, -4, -5, -6, -7, -8, -9, -10, -11, -12][i - 1], amount=1)
            self.operations.append(o)

            o = SymmetryOperation(n=self.n, name="σd", transform_p=lambda i: [2, 1, 6, 5, 4, 3, 8, 7, 12, 11, 10, 9][i - 1], amount=1)
            self.operations.append(o)
            o = SymmetryOperation(n=self.n, name="σd", transform_p=lambda i: [4, 3, 2, 1, 6, 5, 10, 9, 8, 7, 12, 11][i - 1], amount=1)
            self.operations.append(o)
            o = SymmetryOperation(n=self.n, name="σd", transform_p=lambda i: [6, 5, 4, 3, 2, 1, 12, 11, 10, 9, 8, 7][i - 1], amount=1)
            self.operations.append(o)

            o = SymmetryOperation(n=self.n, name="σv", transform_p=lambda i: [1, 6, 5, 4, 3, 2, 7, 12, 11, 10, 9, 8][i - 1], amount=1)
            self.operations.append(o)
            o = SymmetryOperation(n=self.n, name="σv", transform_p=lambda i: [3, 2, 1, 6, 5, 4, 9, 8, 7, 12, 11, 10][i - 1], amount=1)
            self.operations.append(o)
            o = SymmetryOperation(n=self.n, name="σv", transform_p=lambda i: [5, 4, 3, 2, 1, 6, 11, 10, 9, 8, 7, 12][i - 1], amount=1)
            self.operations.append(o)

            self.dipole_operator_symmetry = "A2u"  # 1 1 1 1 -1 -1 -1 -1 -1 1 1, = z
        elif self.n == 6:
            o = SymmetryOperation(n=self.n, name="C6", transform_p=lambda i: [6,1,2,3,4,5][i - 1], amount=1)
            self.operations.append(o)
            o = SymmetryOperation(n=self.n, name="C6", transform_p=lambda i: [2,3,4,5,6,1][i - 1], amount=1)
            self.operations.append(o)

            o = SymmetryOperation(n=self.n, name="C3", transform_p=lambda i: [5,6,1,2,3,4][i - 1], amount=1)
            self.operations.append(o)
            o = SymmetryOperation(n=self.n, name="C3", transform_p=lambda i: [3,4,5,6,1,2][i - 1], amount=1)
            self.operations.append(o)

            o = SymmetryOperation(n=self.n, name="C2", transform_p=lambda i: [4,5,6,1,2,3][i - 1], amount=1)
            self.operations.append(o)

            o = SymmetryOperation(n=self.n, name="C\'2", transform_p=lambda i: [-1,-6,-5,-4,-3,-2][i - 1], amount=1)
            self.operations.append(o)
            o = SymmetryOperation(n=self.n, name="C\'2", transform_p=lambda i: [-3,-2,-1,-6,-5,-4][i - 1], amount=1)
            self.operations.append(o)
            o = SymmetryOperation(n=self.n, name="C\'2", transform_p=lambda i: [-5,-4,-3,-2,-1,-6][i - 1], amount=1)
            self.operations.append(o)

            o = SymmetryOperation(n=self.n, name="C\'\'2", transform_p=lambda i: [-2,-1,-6,-5,-4,-3][i - 1], amount=1)
            self.operations.append(o)
            o = SymmetryOperation(n=self.n, name="C\'\'2", transform_p=lambda i: [-4,-3,-2,-1,-6,-5][i - 1], amount=1)
            self.operations.append(o)
            o = SymmetryOperation(n=self.n, name="C\'\'2", transform_p=lambda i: [-6,-5,-4,-3,-2,-1][i - 1], amount=1)
            self.operations.append(o)

            o = SymmetryOperation(n=self.n, name="i", transform_p=lambda i: [-4,-5,-6,-1,-2,-3][i - 1], amount=1)
            self.operations.append(o)

            o = SymmetryOperation(n=self.n, name="S3", transform_p=lambda i: [-5,-6,-1,-2,-3,-4][i - 1], amount=1)
            self.operations.append(o)
            o = SymmetryOperation(n=self.n, name="S3", transform_p=lambda i: [-3,-4,-5,-6,-1,-2][i - 1], amount=1)
            self.operations.append(o)

            o = SymmetryOperation(n=self.n, name="S6", transform_p=lambda i: [-6,-1,-2,-3,-4,-5][i - 1], amount=1)
            self.operations.append(o)
            o = SymmetryOperation(n=self.n, name="S6", transform_p=lambda i: [-2,-3,-4,-5,-6,-1][i - 1], amount=1)
            self.operations.append(o)

            o = SymmetryOperation(n=self.n, name="σh", transform_p=lambda i: [-1,-2,-3,-4,-5,-6][i - 1], amount=1)
            self.operations.append(o)

            o = SymmetryOperation(n=self.n, name="σd", transform_p=lambda i: [2,1,6,5,4,3][i - 1], amount=1)
            self.operations.append(o)
            o = SymmetryOperation(n=self.n, name="σd", transform_p=lambda i: [4,3,2,1,6,5][i - 1], amount=1)
            self.operations.append(o)
            o = SymmetryOperation(n=self.n, name="σd", transform_p=lambda i: [6,5,4,3,2,1][i - 1], amount=1)
            self.operations.append(o)

            o = SymmetryOperation(n=self.n, name="σv", transform_p=lambda i: [1,6,5,4,3,2][i - 1], amount=1)
            self.operations.append(o)
            o = SymmetryOperation(n=self.n, name="σv", transform_p=lambda i: [3,2,1,6,5,4][i - 1], amount=1)
            self.operations.append(o)
            o = SymmetryOperation(n=self.n, name="σv", transform_p=lambda i: [5,4,3,2,1,6][i - 1], amount=1)
            self.operations.append(o)

            self.dipole_operator_symmetry = "A2u" # 1 1 1 1 -1 -1 -1 -1 -1 1 1, = z
        else:
            raise Exception("unknown number of transforming orbitals")

    def set_up_irreducible_representations(self):
        # if self.n == 6:
            self.irreducible_representations = []
            symmetry_names = list(dict.fromkeys(op.name for op in self.operations))

            i = IrreducibleRepresentation({j: 1 for j in symmetry_names}, name="A1g")
            self.irreducible_representations.append(i)
            A2g_chars = [1, 1, 1, 1, -1, -1 ,1, 1, 1, 1, -1, -1]
            i = IrreducibleRepresentation(dict(zip(symmetry_names, A2g_chars)), name="A2g")
            self.irreducible_representations.append(i)

            B1g_chars = [1, -1, 1, -1, 1, -1, 1,-1, 1, -1, 1, -1 ]
            i = IrreducibleRepresentation(dict(zip(symmetry_names, B1g_chars)), name="B1g")
            self.irreducible_representations.append(i)

            B2g_chars = [1, -1, 1, -1, -1, 1, 1, -1, 1, -1, -1, 1]
            i = IrreducibleRepresentation(dict(zip(symmetry_names, B2g_chars)), name="B2g")
            self.irreducible_representations.append(i)

            E1g_chars = [2, 1, -1, -2, 0, 0, 2, 1, -1, -2, 0, 0]
            i = IrreducibleRepresentation(dict(zip(symmetry_names, E1g_chars)), name="E1g", dimension=2)
            self.irreducible_representations.append(i)

            E2g_chars = [2, -1, -1, 2, 0, 0, 2, -1, -1, 2, 0, 0]
            i = IrreducibleRepresentation(dict(zip(symmetry_names, E2g_chars)), name="E2g", dimension=2)
            self.irreducible_representations.append(i)

            A1u_chars = [1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1, -1]
            i = IrreducibleRepresentation(dict(zip(symmetry_names, A1u_chars)), name="A1u")
            self.irreducible_representations.append(i)

            A2u_chars = [1, 1, 1, 1, -1, -1, -1, -1, -1, -1, 1, 1]
            i = IrreducibleRepresentation(dict(zip(symmetry_names, A2u_chars)), name="A2u")
            self.irreducible_representations.append(i)

            B1u_chars = [1, -1, 1, -1, 1, -1, -1, 1, -1, 1, -1, 1]
            i = IrreducibleRepresentation(dict(zip(symmetry_names, B1u_chars)), name="B1u")
            self.irreducible_representations.append(i)

            B2u_chars = [1, -1, 1, -1, -1, 1, -1, 1, -1, 1, 1, -1]
            i = IrreducibleRepresentation(dict(zip(symmetry_names, B2u_chars)), name="B2u")
            self.irreducible_representations.append(i)

            E1u_chars = [2, 1, -1, -2, 0, 0, -2, -1, 1, 2, 0, 0]
            i = IrreducibleRepresentation(dict(zip(symmetry_names, E1u_chars)), name="E1u", dimension=2)
            self.irreducible_representations.append(i)

            E2u_chars = [2, -1, -1, 2, 0, 0, -2, 1, 1, -2, 0, 0 ]
            i = IrreducibleRepresentation(dict(zip(symmetry_names, E2u_chars)), name="E2u", dimension=2)
            self.irreducible_representations.append(i)
        # else:
        #     raise Exception("unknown number of transforming orbitals")




    def multiply(self, irred_1: str, irred_2: str) -> list[str]:
        if irred_1 == "A1g":
            return [irred_2]
        if irred_2 == "A1g":
            return [irred_1]
        #
        if irred_1 == irred_2:
            if irred_1 == "E1g" or irred_1 == "E2g" or irred_1 == "E1u" or  irred_1 == "E2u":
                return ["A1g", "A2g", "E2g"]
            return ["A1g"]

        if sorted([irred_1, irred_2]) == sorted(["E2g", "E1g"]) or sorted([irred_1, irred_2]) == sorted(["E1u", "E2u"]):
            return ["B1g", "B2g", "E1g"]
        if sorted([irred_1, irred_2]) == sorted(["E1g", "E1u"]) or sorted([irred_1, irred_2]) == sorted(["E2g", "E2u"]):
            return ["A1u", "A2u", "E2u"]
        if sorted([irred_1, irred_2]) == sorted(["E1g", "E2u"]) or sorted([irred_1, irred_2]) == sorted(["E2g", "E1u"]):
            return ["B1u", "B2u", "E1u"]
        # if irred_1 == "E1g":
        #     if irred_2 in ["A2g", "A1g"]:
        #         return ["E1g"]
        #     if irred_2 and ["B1g", "B2g"]:
        #         return ["E2g"]
        #     if irred_2 in ["A1u", "A2u"]:
        #         return ["E1u"]
        #     if  irred_2 in ["B1u", "B2u"]:
        #         return ["E2u"]
        # if irred_2 == "E1g":
        #     if irred_1 in ["A2g", "A1g"]:
        #         return ["E1g"]
        #     if irred_1 and ["B1g", "B2g"]:
        #         return ["E2g"]
        #     if irred_1 in ["A1u", "A2u"]:
        #         return ["E1u"]
        #     if irred_1 in ["B1u", "B2u"]:
        #         return ["E2u"]

        if irred_1 == "E2g":
            if irred_2 in ["A1u", "A2u"]:
                return ["E2u"]
            if irred_2 in ["B1u", "B2u"]:
                return ["E1u"]
        if irred_2 == "E2g":
            if irred_1 in ["A1u", "A2u"]:
                return ["E2u"]
            if irred_1 in ["B1u", "B2u"]:
                return ["E1u"]


        if irred_1 == "A2g":
            if "E" in irred_2:
                return [irred_2]
            return [self.replace_number(irred_2)]
        if irred_2 == "A2g":
            if "E" in irred_1:
                return [irred_1]
            return [self.replace_number(irred_1)]
        if irred_1 == "B1g":
            if irred_2 in ["E1g", "E2g", "E1u", "E2u"]:
                return [self.replace_number(irred_2)]
            return [self.replace_symbol(irred_2)]
        if irred_2 == "B1g":
            if irred_1 in ["E1g", "E2g", "E1u", "E2u"]:
                return [self.replace_number(irred_1)]
            return [self.replace_symbol(irred_1)]
        if irred_1 == "A1u":
            return [self.replace_gu(irred_2)]
        if irred_2 == "A1u":
            return [self.replace_gu(irred_1)]
        if irred_1 == "A2u":
            if irred_2 in ["E1g", "E2g", "E1u", "E2u"]:
                return [self.replace_gu(irred_2)]
            return [self.replace_gu(self.replace_number(irred_2))]
        if irred_2 == "A2u":
            if irred_1 in ["E1g", "E2g", "E1u", "E2u"]:
                return [self.replace_gu(irred_1)]
            return [self.replace_gu(self.replace_number(irred_1))]


        if irred_1 == "B2g":
            if irred_2 in ["E1g", "E2g", "E1u", "E2u"]:
                return [self.replace_number(irred_2)]
            return [self.replace_number(self.replace_symbol(irred_2))]
        if irred_2 == "B2g":
            if irred_1 in ["E1g", "E2g", "E1u", "E2u"]:
                return [self.replace_number(irred_1)]
            return [self.replace_number(self.replace_symbol(irred_1))]



        if irred_1 == "B1u":
            if irred_2 in ["E1g", "E2g", "E1u", "E2u"]:
                return [self.replace_gu(self.replace_number(irred_2))]
            return [self.replace_gu(self.replace_symbol(irred_2))]
        if irred_2 == "B1u":
            if irred_1 in ["E1g", "E2g", "E1u", "E2u"]:
                return [self.replace_gu(self.replace_number(irred_1))]
            return [self.replace_gu(self.replace_symbol(irred_1))]
        if irred_1 == "B2u":
            if irred_2 in ["E1g", "E2g", "E1u", "E2u"]:
                return [self.replace_gu(self.replace_number(irred_2))]
            return [self.replace_gu(self.replace_symbol(self.replace_number(irred_2)))]
        if irred_2 == "B2u":
            if irred_1 in ["E1g", "E2g", "E1u", "E2u"]:
                return [self.replace_gu(self.replace_number(irred_1))]
            return [self.replace_gu(self.replace_symbol(self.replace_number(irred_1)))]



if __name__ == '__main__':
    # D6h().print_character_table()

    D6h().print_multiplication_table()