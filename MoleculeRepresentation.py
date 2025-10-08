import itertools
from fractions import Fraction
import sympy as sp


from IrreducibleRepresentation import IrreducibleRepresentation
from PointGroups.D4h import D4h
from SALC import SALC, norm_and_group_SALCs
from SymmetryOperation import SymmetryOperation
from is_multiple import is_multiple
from tst.solve_saekular_equation import calculate_hueckel_secular_equation


class MoleculeRepresentation():
    def __init__(self, n:int):
        self.n = n
        self.circular = True
        self.set_up_point_group()


    def set_up_point_group(self):
        if self.n == 4 and self.circular:# D4h, cyclopentadiene
            self.pointgroup = D4h()
        elif self.n == 6:
            pass
        else:
            raise Exception("Not implemented")

    def set_to_s_orbitals(self):
        # can't be undone!!!
        self.s_orbital_active = True
        for s in range(len(self.operations)):
            old_transform = self.operations[s].transform_p
            self.operations[s].transform_p = lambda x, f=old_transform: abs(f(x))

    def set_up_irreducible_representations(self):
        symmetry_names = list(dict.fromkeys(op.name for op in self.operations))
        if self.n == 4:
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
        elif self.n == 6:
            print("TODO")
        else:
            raise Exception("Not implemented")


    def get_reducible_representation_for_ring_p_orbitals(self):
        reducible_representation = {}
        for op in self.pointgroup.operations:
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

    # def get_symmetry_operation_by_name(self, name:str):#INFO Impossible, because names can be duplicate

    def get_irreducible_representation_by_name(self, name:str):
        for i in self.pointgroup.irreducible_representations:
            if i.name == name:
                return i
        raise Exception("Not found")

    def decomposing_into_irreducible_representations(self, reducible_representation:dict):
        """
        Zerlegt die reduzible Darstellung in irreduzible Darstellungen.
        a_i =  (1/n))  * sum(   Xi_red(R) * Xi_irred^(i) (R)^*   )
        """
        decomposition = {}

        for irrep in self.pointgroup.irreducible_representations:
            a_i = 0
            for op in self.pointgroup.operations:
                #summiert nicht mehrfach über analoge sym-operationen, daher zusätzlicher Faktor in Berechnung
                a_i += reducible_representation[op.name] * irrep.characters[op.name] * op.amount
            a_i /= self.pointgroup.group_order()
            assert a_i - int(a_i) < 1e-10
            decomposition[irrep.name] = int(a_i)
        return decomposition


    def get_expected_number_of_SALC_per_irreducible_representation(self, irreducible_representation:dict ):
        number_of_SALC_per_irred = {}
        for i in self.pointgroup.irreducible_representations:
            number_of_SALC_per_irred[i.name] = i.dimension * irreducible_representation[i.name]
        return number_of_SALC_per_irred

    def project(self, irreducible_representation:dict, p_orbital_index:int):
        # Symbolische Basisfunktionen
        p_symbols = sp.symbols(f"p1:{self.n + 1}")
        SALCs = []
        for irred_name, value in irreducible_representation.items():
            eq_expr = 0
            collected_p_orbitals_in_expression = {}
            if value != 0:
                irred = self.get_irreducible_representation_by_name(irred_name)
                for op in self.pointgroup.operations:
                    op_name = op.name
                    op_value = irred.characters[op_name] if op_name in irred.characters.keys() else 0
                # for op_name, op_value in irred.characters.items():
                    if op_value != 0:
                        # op = self.get_symmetry_operation_by_name(op_name)
                        transformed_p = op.transform_p(p_orbital_index)
                        # print(p_orbital_index, "->", transformed_p)
                        coeff = op_value if transformed_p > 0 else -op_value
                        transformed_p = p_symbols[abs(transformed_p) - 1]
                        eq_expr += op.amount * coeff * transformed_p
                        collected_p_orbitals_in_expression[transformed_p] = collected_p_orbitals_in_expression.get(transformed_p, 0) + op.amount * coeff
                        # print("adding", op.name, "=",coeff , f"{transformed_p}")
                eq_expr /= self.pointgroup.group_order()
                collected_p_orbitals_in_expression = {k: Fraction(v, self.pointgroup.group_order()) for k, v in collected_p_orbitals_in_expression.items()}
                s = SALC(n=self.n, irred=irred_name, p_orbital_prefactors=collected_p_orbitals_in_expression, equation=eq_expr)
                # s.print()
                SALCs.append(s)
        return SALCs



    def get_all_SALCs(self):
        reducible_representation = self.get_reducible_representation_for_ring_p_orbitals()
        # print(reducible_representation)

        irreducible_representation = self.decomposing_into_irreducible_representations(reducible_representation)
        # print(irreducible_representation)

        all_SALCs = []
        for n in range(1, self.n + 1):
            SALCs = self.project(irreducible_representation, p_orbital_index=n)
            for s in SALCs:
                # s.print()
                checking_list = [is_multiple(x.equation, s.equation) for x in all_SALCs]# negative version of LC is still the same linear combination
                # old_list = [old_s.prefactors_of_AOs for old_s in all_SALCs]
                if not True in checking_list:
                    # print("not in", s.prefactors_of_AOs, s.irred)
                    all_SALCs.append(s)

        # print(len(all_SALCs), irreducible_representation, sum(irreducible_representation.values()))

        # expectation:
        assert len(all_SALCs) == sum(self.get_expected_number_of_SALC_per_irreducible_representation(irreducible_representation).values())
        return all_SALCs



    def orbitals_adjoint(self, p1, p2):
        alpha, beta = sp.symbols("alpha beta")
        if p1 == p2:
            return alpha
        if p1 + 1 == p2 or p1 -1 == p2 :
            return beta
        if self.circular and (
                            (p1-1 == 0 and p2 == self.n) or
                            (p2-1 == 0 and p1 == self.n)  ):
            return beta
        return 0

    def h_eff(self, salc_1: SALC, salc_2: SALC):
        if not salc_1.normalized or not salc_2.normalized:
            raise Exception("not normalized salcs!")
        h_eff_integral = 0
        for p_orbital in range(1, len(salc_1.prefactors_of_AOs) + 1):
            if salc_1.prefactors_of_AOs[p_orbital - 1] != 0:
                part = 0
                for combi in range(1, len(salc_2.prefactors_of_AOs) + 1):
                    if salc_2.prefactors_of_AOs[combi - 1] != 0:
                        factor = salc_1.prefactors_of_AOs[p_orbital - 1] * salc_2.prefactors_of_AOs[combi - 1]
                        part += factor * self.orbitals_adjoint(p_orbital, combi)
                        print("p", p_orbital, "|", "p", combi, "=")
                        print("+", self.orbitals_adjoint(p_orbital, combi) ,"*", factor)
                h_eff_integral += part
        return h_eff_integral

    def get_effective_hamilton_matrix(self, SALCs: list[SALC], irred:str):
        # print("get_effective_hamilton_matrix", [s.irred for s in SALCs], flush=True)
        # print(f"\n=== {irred} ===", type(irred), flush=True)

        n = len(SALCs)
        H = sp.zeros(n)
        H_show = sp.zeros(n)
        if len(SALCs) == 1:
                h = sp.Symbol(f"H_{irred}")
                # print("H =", h)
                H[0, 0] = self.h_eff(salc_1=SALCs[0], salc_2=SALCs[0])
        else:
                for i, j in itertools.product(range(n), repeat=2):
                    H[i, j] = self.h_eff(salc_1=SALCs[i], salc_2=SALCs[j])
                    H_show[i, j] = sp.Symbol(f"H_{irred}_{i}{j}")
                # print("H =")
                # sp.pprint(H_show)
                # print("=")
                # sp.pprint(H)
        return H

    def get_energy_levels(self):
        SALCs = self.get_all_SALCs()
        SALCs_by_irred = norm_and_group_SALCs(SALCs)

        result = []
        alpha, beta = sp.symbols(f"alpha beta")
        sorting_dict_values = {alpha: 0, beta: -1}
        for irred, salcs in SALCs_by_irred.items():
            H = self.get_effective_hamilton_matrix(SALCs=salcs, irred=irred)
            result_irred = calculate_hueckel_secular_equation(H, info=irred, sorting_dict_values=sorting_dict_values)
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
    p = MoleculeRepresentation(n=4)
    p.set_to_s_orbitals()

    reducible_representation = p.get_reducible_representation_for_ring_p_orbitals()
    print(reducible_representation)

    irreducible_representation = p.decomposing_into_irreducible_representations(reducible_representation)
    print(irreducible_representation)

    salcs = p.get_all_SALCs()
    for s in salcs:
        print(s.irred, s.equation)

    # p.project(irreducible_representation, p_orbital_index=1)
    p.set_to_s_orbitals()
    s_mo_orbitals = p.get_energy_levels()

    for s in s_mo_orbitals:
        print("\t", s.symmetry, ":\t", s.eigenvalue)
