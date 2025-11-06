from IrreducibleRepresentation import IrreducibleRepresentation
from PointGroups.PointGroup import PointGroup
from SymmetryOperation import SymmetryOperation


class D2h(PointGroup):

    def __init__(self, n:int=4):
        super().__init__(order=8, n=n)
        self.total_symmetric_representation = "Ag"
        self.set_up_symmetry_operations()
        self.set_up_irreducible_representations()

    def set_up_symmetry_operations(self):
        o = SymmetryOperation(n=self.n, name="E", transform_p=lambda i: i, amount = 1)
        self.operations.append(o)

        o = SymmetryOperation(n=self.n, name="C2(z)", transform_p=lambda i: [][i - 1], amount=1)
        self.operations.append(o)
        o = SymmetryOperation(n=self.n, name="C2(y)", transform_p=lambda i: [][i - 1], amount=1)
        self.operations.append(o)
        o = SymmetryOperation(n=self.n, name="C2(x)", transform_p=lambda i: [][i - 1], amount=1)
        self.operations.append(o)

        o = SymmetryOperation(n=self.n, name="i", transform_p=lambda i: [][i - 1], amount=1)
        self.operations.append(o)

        o = SymmetryOperation(n=self.n, name="σ(xy)", transform_p=lambda i: [][i - 1], amount=1)
        self.operations.append(o)
        o = SymmetryOperation(n=self.n, name="σ(xz)", transform_p=lambda i: [][i - 1], amount=1)
        self.operations.append(o)
        o = SymmetryOperation(n=self.n, name="σ(yz)", transform_p=lambda i: [][i - 1], amount=1)
        self.operations.append(o)


        self.dipole_operator_symmetry = ""


    def set_up_irreducible_representations(self):
        symmetry_names = list(dict.fromkeys(op.name for op in self.operations))
        i = IrreducibleRepresentation(characters={j: 1 for j in symmetry_names}, name="Ag")
        self.irreducible_representations.append(i)

        B1g_chars = [1, 1, -1, -1, 1, 1, -1, -1]
        i = IrreducibleRepresentation(characters=dict(zip(symmetry_names, B1g_chars)), name="B1g")
        self.irreducible_representations.append(i)

        B2g_chars = [1,	-1,	1,	-1,	1,	-1,	1,	-1]
        i = IrreducibleRepresentation(characters=dict(zip(symmetry_names, B2g_chars)), name="B2g")
        self.irreducible_representations.append(i)

        B3g_chars = [1,	-1,	-1,	1,	1,	-1,	-1,	1]
        i = IrreducibleRepresentation(characters=dict(zip(symmetry_names, B3g_chars)), name="B3g")
        self.irreducible_representations.append(i)
        Au_chars = [1,	1,	1,	1,	-1,	-1,	-1,	-1]
        i = IrreducibleRepresentation(characters=dict(zip(symmetry_names, Au_chars)), name="Au")
        self.irreducible_representations.append(i)
        B1u_chars = [1,	1,	-1,	-1,	-1,	-1,	1,	1]
        i = IrreducibleRepresentation(characters=dict(zip(symmetry_names, B1u_chars)), name="B1u")
        self.irreducible_representations.append(i)
        B2u_chars = [1,	-1,	1,	-1,	-1,	1,	-1,	1]
        i = IrreducibleRepresentation(characters=dict(zip(symmetry_names, B2u_chars)), name="B2u")
        self.irreducible_representations.append(i)

        B3u_chars = [1,	-1,	-1,	1,	-1,	1,	1,	-1]
        i = IrreducibleRepresentation(characters=dict(zip(symmetry_names, B3u_chars)), name="B3u")
        self.irreducible_representations.append(i)


    def multiply(self, irred_1:str, irred_2:str) -> list[str]:
        possible_irreducible_representations = [i.name for i in self.irreducible_representations]
        if irred_1 not in possible_irreducible_representations or irred_2 not in possible_irreducible_representations:
            raise Exception(f"unknown irreducible representation {irred_1} or {irred_2}")

        if irred_1 == "Ag":
            return [irred_2]
        if irred_2 == "Ag":
            return [irred_1]
        if irred_1 == irred_2:
            return ["Ag"]

        if irred_1 == "Au":
            return [self.replace_gu(irred_2)]
        if irred_2 == "Au":
            return [self.replace_gu(irred_1)]

        if irred_2.replace("g","u") == irred_1:
            return ["Au"]


        if ("g" in irred_1 and "g" in irred_2) or ("u" in irred_1 and "u" in irred_2):
            gu = "g"
        else:
            gu = "u"
        if ("1" in irred_1 and "1" in irred_2) or ("2" in irred_1 and "2" in irred_2) or ("3" in irred_1 and "3" in irred_2):
            return ["A"+gu]
        if ("1" in irred_1 and "2" in irred_2) or ("2" in irred_1 and "1" in irred_2):
            return ["B3"+gu]
        if ("1" in irred_1 and "3" in irred_2) or ("3" in irred_1 and "1" in irred_2):
            return ["B2"+gu]
        if ("3" in irred_1 and "2" in irred_2) or ("2" in irred_1 and "3" in irred_2):
            return ["B1"+gu]



        # if irred_1 == "B1g":
        #     if irred_2 == "B2g":
        #         return ["B3g"]
        #     if irred_2 == "B3g":
        #         return ["B2g"]
        #     # if irred_2 == "B1u":
        #     #     return ["Au"]
        #     if irred_2 == "B2u":
        #         return ["B3u"]
        #     if irred_2 == "B3u":
        #         return ["B2u"]
        # if irred_2 == "B1g":
        #     if irred_1 == "B2g":
        #         return ["B3g"]
        #     if irred_1 == "B3g":
        #         return ["B2g"]
        #     if irred_1 == "B2u":
        #         return ["B3u"]
        #     if irred_1 == "B3u":
        #         return ["B2u"]
        #
        # if irred_1 == "B2g":
        #     if irred_2 == "B3g":
        #         return ["B1g"]
        #     if irred_2 == "B1u":
        #         return ["B3u"]
        #     if irred_2 == "B3u":
        #         return ["B1u"]




        raise Exception("nothing should remain")




if __name__ == "__main__":
    D2h = D2h()
    # D2h.print_character_table()
    D2h.print_multiplication_table()