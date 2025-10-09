from IrreducibleRepresentation import IrreducibleRepresentation
from PointGroups.PointGroup import PointGroup
from SymmetryOperation import SymmetryOperation


class D4h(PointGroup):

    def __init__(self, n:int=4):
        super().__init__(order=16, n=n)
        self.total_symmetric_representation = "A1g"
        self.set_up_symmetry_operations()
        self.set_up_irreducible_representations()

    def set_up_symmetry_operations(self):
        o = SymmetryOperation(n=self.n, name="E", transform_p=lambda i: i, amount = 1)
        self.operations.append(o)
        o = SymmetryOperation(n=self.n, name="C4(z)", transform_p=lambda i: [4,1,2,3][i-1], amount=1)

        self.operations.append(o)
        o = SymmetryOperation(n=self.n, name="C4(z)", transform_p=lambda i: [2, 3, 4, 1][i - 1], amount=1)
        self.operations.append(o)

        o = SymmetryOperation(n=self.n, name="C2", transform_p=lambda i: [3, 4, 1, 2][i - 1], amount=1)
        self.operations.append(o)

        o = SymmetryOperation(n=self.n, name="C\'2", transform_p=lambda i: [-1, -4, -3, -2][i - 1], amount=1)
        self.operations.append(o)
        o = SymmetryOperation(n=self.n, name="C\'2", transform_p=lambda i: [-3, -2, -1, -4][i - 1], amount=1)
        self.operations.append(o)

        o = SymmetryOperation(n=self.n, name="2C\'\'2", transform_p=lambda i: [-2, -1, -4, -3][i - 1], amount=1)
        self.operations.append(o)
        o = SymmetryOperation(n=self.n, name="2C\'\'2", transform_p=lambda i: [-4, -3, -2, -1][i - 1], amount=1)
        self.operations.append(o)

        o = SymmetryOperation(n=self.n, name="i", transform_p=lambda i: [-3, -4, -1, -2][i - 1], amount=1)
        self.operations.append(o)

        o = SymmetryOperation(n=self.n, name="S4", transform_p=lambda i: [-4, -1, -2, -3][i - 1], amount=1)
        self.operations.append(o)
        o = SymmetryOperation(n=self.n, name="S4", transform_p=lambda i: [-2, -3, -4, -1][i - 1], amount=1)
        self.operations.append(o)

        o = SymmetryOperation(n=self.n, name="σh", transform_p=lambda i: [-1, -2, -3, -4][i - 1], amount=1)
        self.operations.append(o)

        o = SymmetryOperation(n=self.n, name="σv", transform_p=lambda i: [1, 4, 3, 2][i - 1], amount=1)
        self.operations.append(o)
        o = SymmetryOperation(n=self.n, name="σv", transform_p=lambda i: [3, 2, 1, 4][i - 1], amount=1)
        self.operations.append(o)

        o = SymmetryOperation(n=self.n, name="σd", transform_p=lambda i: [2, 1, 4, 3][i - 1], amount=1)
        self.operations.append(o)
        o = SymmetryOperation(n=self.n, name="σd", transform_p=lambda i: [4, 3, 2, 1][i - 1], amount=1)
        self.operations.append(o)

        self.dipole_operator_symmetry = "A2u"


    def set_up_irreducible_representations(self):
        symmetry_names = list(dict.fromkeys(op.name for op in self.operations))
        i = IrreducibleRepresentation(characters={j: 1 for j in symmetry_names}, name="A1g")
        self.irreducible_representations.append(i)
        A2g_chars = [1	,1,	1,	-1	,-1	,1,	1,	1,	-1,	-1	]
        i = IrreducibleRepresentation(characters=dict(zip(symmetry_names, A2g_chars)), name="A2g")
        self.irreducible_representations.append(i)
        B1g_chars = [1,	-1,	1,	1,	-1	,1	,-1	,1,	1,	-1]
        i = IrreducibleRepresentation(characters=dict(zip(symmetry_names, B1g_chars)), name="B1g")
        self.irreducible_representations.append(i)
        B2g_chars = [1	,-1	,1	,-1	,1,	1,	-1,	1,	-1,	1	]
        i = IrreducibleRepresentation(characters=dict(zip(symmetry_names, B2g_chars)), name="B2g")
        self.irreducible_representations.append(i)
        Eg_chars = [2,	0,	-2,	0,	0,	2,	0,	-2,	0,	0]
        i = IrreducibleRepresentation(characters=dict(zip(symmetry_names, Eg_chars)), name="Eg", dimension=2)
        self.irreducible_representations.append(i)
        A1u_chars = [1,	1,	1,	1,	1,	-1,	-1,	-1,	-1,	-1	]
        i = IrreducibleRepresentation(characters=dict(zip(symmetry_names, A1u_chars)), name="A1u")
        self.irreducible_representations.append(i)
        A2u_chars = [1,	1	,1,	-1,	-1,	-1	,-1,	-1	,1,	1]
        i = IrreducibleRepresentation(characters=dict(zip(symmetry_names, A2u_chars)), name="A2u")
        self.irreducible_representations.append(i)
        B1u_chars = [1,	-1,	1,	1,	-1	,-1,	1	,-1	,-1	,1	]
        i = IrreducibleRepresentation(characters=dict(zip(symmetry_names, B1u_chars)), name="B1u")
        self.irreducible_representations.append(i)
        B2u_chars = [1,	-1	,1	,-1,	1,	-1,	1	,-1	,1,	-1	]
        i = IrreducibleRepresentation(characters=dict(zip(symmetry_names, B2u_chars)), name="B2u")
        self.irreducible_representations.append(i)
        Eu_chars = [2,	0, -2	,0	,0,	-2	,0,	2,	0,	0]
        i = IrreducibleRepresentation(characters=dict(zip(symmetry_names, Eu_chars)), name="Eu", dimension=2)
        self.irreducible_representations.append(i)


    def multiply(self, irred_1:str, irred_2:str) -> list[str]:
        possible_irreducible_representations = [i.name for i in self.irreducible_representations]
        if irred_1 not in possible_irreducible_representations or irred_2 not in possible_irreducible_representations:
            raise Exception(f"unknown irreducible representation {irred_1} or {irred_2}")

        if irred_1 == "A1g":
            return [irred_2]
        if irred_2 == "A1g":
            return [irred_1]
        if irred_1 == irred_2:
            if irred_1 == "Eg" or irred_1 == "Eu":
                return ["A1g", "A2g", "B1g", "B2g"]
            else:
                return ["A1g"]
        if (irred_1 == "Eu" and irred_2 == "Eg") or (irred_1 == "Eg" and irred_2 == "Eu"):
            return ["A1u", "A2u", "B1u", "B2u"]
        if irred_1 == "Eu":
            return ["Eg"] if "u" in irred_2 else ["Eu"]
        if irred_2 == "Eu":
            return ["Eg"] if "u" in irred_1 else ["Eu"]

        if irred_1 == "Eg":
            return ["Eg"] if "g" in irred_2 else ["Eu"]
        if irred_2 == "Eg":
            return ["Eg"] if "g" in irred_1 else ["Eu"]

        if irred_1 == "A2g":
            return [irred_2.replace("2","z").replace("1","2").replace("z","1")]
        if irred_2 == "A2g":
            return [irred_1.replace("2","z").replace("1","2").replace("z","1")]

        if irred_1 == "B1g":
            return [irred_2.replace("A","z").replace("B","A").replace("z","B")]
        if irred_2 == "B1g":
            return [irred_1.replace("A","z").replace("B","A").replace("z","B")]

        if irred_1 == "A1u":
            return [irred_2.replace("g","z").replace("u","g").replace("z","u")]
        if irred_2 == "A1u":
            return [irred_1.replace("g","z").replace("u","g").replace("z","u")]

        if irred_1 == "B2g":
            version_1 = irred_2.replace("2","z").replace("1","2").replace("z","1")
            version_2 = version_1.replace("A","z").replace("B","A").replace("z","B")
            return [version_2]
        if irred_2 == "B2g":
            version_1 = irred_1.replace("2","z").replace("1","2").replace("z","1")
            version_2 = version_1.replace("A","z").replace("B","A").replace("z","B")
            return [version_2]

        if irred_1 == "A2u":
            version_1 = irred_2.replace("2","z").replace("1","2").replace("z","1")
            version_2 = version_1.replace("g","z").replace("u","g").replace("z","u")
            return [version_2]
        if irred_2 == "A2u":
            version_1 = irred_1.replace("2","z").replace("1","2").replace("z","1")
            version_2 = version_1.replace("g","z").replace("u","g").replace("z","u")
            return [version_2]
        if irred_1 == "B2u" and irred_2 == "B2u":
            return ["A1g"]
        if (irred_1 == "B2u" and irred_2 == "B1u") or (irred_2 == "B2u" and irred_1 == "B1u"):
            return ["A2g"]
        raise Exception("nothing should remain")




if __name__ == "__main__":
    D4h = D4h()

    D4h.print_multiplication_table()