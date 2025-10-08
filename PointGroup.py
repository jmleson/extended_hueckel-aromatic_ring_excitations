import itertools
from fractions import Fraction
import sympy as sp


from IrreducibleRepresentation import IrreducibleRepresentation
from SALC import SALC, norm_and_group_SALCs
from SymmetryOperation import SymmetryOperation
from tst.solve_saekular_equation import calculate


class PointGroup():
    def __init__(self, n:int):
        self.n = n
        self.operations = []
        self.set_up_symmetry_operations()
        self.irreducible_representations = []
        self.set_up_irreducible_representations()

        self.circular = True

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
            i = IrreducibleRepresentation(characters={j: 1 for j in symmetry_names}, name="A\'1")
            self.irreducible_representations.append(i)
            A2_chars = [1,1,-1,1,1,-1]
            i = IrreducibleRepresentation(characters=dict(zip(symmetry_names, A2_chars)), name="A\'2")
            self.irreducible_representations.append(i)
            E_chars = [2,-1,0,2,-1,0]
            i = IrreducibleRepresentation(characters=dict(zip(symmetry_names, E_chars)), name="E\'")
            self.irreducible_representations.append(i)
            A1_chars = [1,1,1,-1,-1,-1]
            i = IrreducibleRepresentation(characters=dict(zip(symmetry_names, A1_chars)), name="A\'\'1")
            self.irreducible_representations.append(i)
            A2_chars = [1,1,-1,-1,-1,1]
            i = IrreducibleRepresentation(characters=dict(zip(symmetry_names, A2_chars)), name="A\'\'2")
            self.irreducible_representations.append(i)
            E_chars = [2, -1, 0, -2, 1, 0]
            i = IrreducibleRepresentation(characters=dict(zip(symmetry_names, E_chars)), name="E\'\'")
            self.irreducible_representations.append(i)
        elif self.n == 6:
            print("TODO")
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
                return i
        raise Exception("Not found")

    def get_irreducible_representation_by_name(self, name:str):
        for i in self.irreducible_representations:
            if i.name == name:
                return i
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

    def project(self, irreducible_representation:dict, p_orbital_index:int):
        # Symbolische Basisfunktionen
        p_symbols = sp.symbols(f"p1:{self.n + 1}")
        SALCs = []

        for irred_name, value in irreducible_representation.items():
            eq_expr = 0
            collected_p_orbitals_in_expression = {}
            if value != 0:
                irred = self.get_irreducible_representation_by_name(irred_name)

                for op_name, op_value in irred.characters.items():
                    if op_value != 0:
                        op = self.get_symmetry_operation_by_name(op_name)
                        transformed_p = op.transform_p(p_orbital_index)
                        coeff = op_value if transformed_p > 0 else -op_value
                        transformed_p = p_symbols[abs(transformed_p) - 1]
                        eq_expr += coeff * transformed_p
                        collected_p_orbitals_in_expression[transformed_p] = collected_p_orbitals_in_expression.get(transformed_p, 0) + coeff
                eq_expr /= self.group_order()
                collected_p_orbitals_in_expression = {k: Fraction(v, self.group_order()) for k, v in collected_p_orbitals_in_expression.items()}
                s = SALC(n=self.n, irred=irred_name, p_orbital_prefactors=collected_p_orbitals_in_expression, equation=eq_expr)
                # s.print()
                SALCs.append(s)
        return SALCs



    def get_all_SALCs(self):
        reducible_representation = self.get_reducible_representation_for_ring_p_orbitals()
        print(reducible_representation)

        irreducible_representation = self.decomposing_into_irreducible_representations(reducible_representation)
        print(irreducible_representation)

        all_SALCs = []
        for n in range(1, self.n + 1):
            SALCs = self.project(irreducible_representation, p_orbital_index=n)
            for s in SALCs:
                # s.print()
                usual_list = s.prefactors_of_AOs
                negative_list = [-x for x in usual_list]# negative version of LC is still the same linear combination
                old_list = [old_s.prefactors_of_AOs for old_s in all_SALCs]
                if usual_list not in old_list and negative_list not in old_list :
                    # print("not in", s.prefactors_of_AOs)
                    all_SALCs.append(s)

        # print(all_SALCs, len(all_SALCs), irreducible_representation, sum(irreducible_representation.values()))
        assert len(all_SALCs) == sum(irreducible_representation.values())
        return all_SALCs


    def orbitals_adjoint(self, p1, p2):
        alpha, beta = sp.symbols("alpha beta")
        if p1 == p2:
            return alpha
        if p1 + 1 == p2 or p1 -1 == p2 :
            return beta
        if self.circular and p1 -1 == 0 and p2 == self.n:
            return beta
        return 0

    def h_eff(self, salc_1: SALC, salc_2: SALC):
        h_eff_integral = 0
        for p_orbital in range(1, len(salc_1.prefactors_of_AOs) + 1):
            if salc_1.prefactors_of_AOs[p_orbital - 1] != 0:
                part = 0
                for combi in range(1, len(salc_2.prefactors_of_AOs) + 1):
                    if salc_2.prefactors_of_AOs[combi - 1] != 0:
                        factor = salc_1.prefactors_of_AOs[p_orbital - 1] * salc_2.prefactors_of_AOs[combi - 1]
                        part += factor * self.orbitals_adjoint(p_orbital, combi)
                        # print("p", p_orbital, "|", "p", combi, "=")
                        # print("+", p.orbitals_adjoint(p_orbital, combi) ,"*", factor)
                h_eff_integral += part
        return h_eff_integral

    def get_effective_hamilton_matrix(self, SALCs: list[SALC], irred:str):
        print("get_effective_hamilton_matrix", [s.irred for s in SALCs], flush=True)
        print(f"\n=== {irred} ===", type(irred), flush=True)

        if len(SALCs) == 1:
                h = sp.Symbol(f"H_{irred}")
                print("H =", h)
                return [h]
        else:
                n = len(SALCs)
                H = sp.zeros(n)
                H_show = sp.zeros(n)
                for i, j in itertools.product(range(n), repeat=2):
                    H[i, j] = self.h_eff(salc_1=SALCs[i], salc_2=SALCs[j])
                    H_show[i, j] = sp.Symbol(f"H_{irred}_{i}{j}")
                print("H =")
                sp.pprint(H_show)
                print("=")
                sp.pprint(H)
                return H

    def get_energy_levels(self):
        SALCs = self.get_all_SALCs()
        SALCs_by_irred = norm_and_group_SALCs(SALCs)

        result = []
        alpha, beta = sp.symbols(f"alpha beta")
        sorting_dict_values = {alpha: 0, beta: -1}
        for irred, salcs in SALCs_by_irred.items():
            H = self.get_effective_hamilton_matrix(SALCs=salcs, irred=irred)
            result_irred = calculate(H, info=irred, sorting_dict_values=sorting_dict_values)
            for i in result_irred:
                i.symmetry = irred
                result.append(i)
        # sorting:
        molecule_orbitals = sorted(
            result,
            key=lambda m: m.eigenvalue.subs(sorting_dict_values)
        )
        return molecule_orbitals


if __name__ == "__main__":
    p = PointGroup(n=4)
    reducible_representation = p.get_reducible_representation_for_ring_p_orbitals()
    print(reducible_representation)

    irreducible_representation = p.decomposing_into_irreducible_representations(reducible_representation)
    print(irreducible_representation)

    p.project(irreducible_representation, p_orbital_index=1)

