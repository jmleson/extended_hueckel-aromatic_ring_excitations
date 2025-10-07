from fractions import Fraction

from IrreducibleRepresentation import IrreducibleRepresentation
from SymmetryOperation import SymmetryOperation


class PointGroup():
    def __init__(self, n:int):
        self.n = n
        self.operations = []
        self.set_up_symmetry_operations()
        self.irreducible_representations = []
        self.set_up_irreducible_representations()

    def set_up_symmetry_operations(self):
        o = SymmetryOperation(n=self.n, name="E", transform_p=lambda i: i, amount = 1)
        self.operations.append(o)
        if self.n == 3:
            #TODO check
            o = SymmetryOperation(n=self.n, name="C3", transform_p=lambda i: [3, 1, 2][i-1], amount=2)
            self.operations.append(o)
            o = SymmetryOperation(n=self.n, name="C\'2", transform_p=lambda i: [-1,-3,-2][i-1], amount=3)
            self.operations.append(o)
            o = SymmetryOperation(n=self.n, name="σh", transform_p=lambda i: [-1,-2,-3][i-1], amount=1)
            self.operations.append(o)
            o = SymmetryOperation(n=self.n, name = "S3", transform_p=lambda i: [-3,-1,-2][i-1], amount=2)
            self.operations.append(o)
            o = SymmetryOperation(n=self.n, name="σv", transform_p=lambda i: [1, 3, 2][i-1], amount=3)
            self.operations.append(o)
        elif self.n == 6:
            pass
        else:
            raise Exception("Not implemented")

    def set_up_irreducible_representations(self):
        symmetry_names = [op.name for op in self.operations]
        if self.n == 3:
            i = IrreducibleRepresentation({j: 1 for j in symmetry_names}, name="A\'1")
            self.irreducible_representations.append(i)
            A2_chars = [1,1,-1,1,1,-1]
            i = IrreducibleRepresentation(dict(zip(symmetry_names, A2_chars)), name="A\'2")
            self.irreducible_representations.append(i)
            E_chars = [2,-1,0,2,-1,0]
            i = IrreducibleRepresentation(dict(zip(symmetry_names, E_chars)), name="E\'")
            self.irreducible_representations.append(i)
            A1_chars = [1,1,1,-1,-1,-1]
            i = IrreducibleRepresentation(dict(zip(symmetry_names, A1_chars)), name="A\'\'1")
            self.irreducible_representations.append(i)
            A2_chars = [1,1,-1,-1,-1,1]
            i = IrreducibleRepresentation(dict(zip(symmetry_names, A2_chars)), name="A\'\'2")
            self.irreducible_representations.append(i)
            E_chars = [2, -1, 0, -2, 1, 0]
            i = IrreducibleRepresentation(dict(zip(symmetry_names, E_chars)), name="E\'\'")
            self.irreducible_representations.append(i)
        else:
            raise Exception("Not implemented")




    def group_order(self):
        return sum([op.amount for op in self.operations])

    def get_reducible_representation_for_ring_p_orbitals(self):
        reducible_representation = {}
        for op in self.operations:
            sum = 0
            for orbital in range(1, self.n+1):
                transformed_orbital_at_place_orbital = op.transform_p(orbital)
                if transformed_orbital_at_place_orbital == orbital:
                    sum += 1
                elif -transformed_orbital_at_place_orbital == orbital:
                    sum -= 1
            # print(sum, op.name)
            reducible_representation[op.name] = sum
        return reducible_representation

    def get_symmetry_operation_by_name(self, name:str):
        for i in self.operations:
            if i.name == name:
                return name
        raise Exception("Not found")

    def decomposing_into_irreducible_representations(self, reducible_representation:dict):
        """
        Zerlegt die reduzible Darstellung in irreduzible Darstellungen.
        a_i =  (1/n))  * sum(   Xi_red(R) * Xi_irred^(i) (R)^*   )
        """
        decomposition = {}

        for irrep in self.irreducible_representations:
            a_i = 0
            for op in self.operations:
                #summiert nicht mehrfach über analoge sym-operationen, daher zusätzlicher Faktor in Berechnung
                a_i += reducible_representation[op.name] * irrep.characters[op.name] * op.amount
            a_i /= self.group_order()
            assert a_i - int(a_i) < 1e-10
            decomposition[irrep.name] = int(a_i)
        return decomposition


if __name__ == "__main__":
    p = PointGroup(n=2)
    reducible_representation = p.get_reducible_representation_for_ring_p_orbitals()
    print(reducible_representation)

    irred = p.decomposing_into_irreducible_representations(reducible_representation)
    print(irred)
