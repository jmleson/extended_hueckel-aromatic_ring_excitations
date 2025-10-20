import itertools
from collections import defaultdict
from fractions import Fraction
import sympy as sp

from PointGroups.C2v import C2v
from PointGroups.D4h import D4h
from PointGroups.D6h import D6h
from SALC import SALC, norm_and_group_SALCs, linear_independent
from is_multiple import is_multiple
from tst.ortho import get_linear_independent_SALCs
from tst.solve_saekular_equation import calculate_hueckel_secular_equation


class MoleculeRepresentation():
    def __init__(self, n:int, bound_cl_to_c_positions:list[int]):
        if not (len(bound_cl_to_c_positions) == 0 or
                (len(bound_cl_to_c_positions) == 1 and bound_cl_to_c_positions[0] == 1)):
            raise Exception("possibility not yet implemented")
        self.n = n
        self.circular = True
        self.bound_cl_to_c_positions = bound_cl_to_c_positions
        self.s_orbital_active = False
        self.set_up_point_group()


    def set_up_point_group(self):
        if len(self.bound_cl_to_c_positions) == 0:
            if self.n == 4 and not self.circular:
                # bend cyclopentadiene
                self.point_group = C2v()
            if self.n == 4 and self.circular:# D4h, cyclopentadiene
                self.pointgroup = D4h()
            elif self.n == 6 and self.circular:
                self.pointgroup = D6h()
            else:
                raise Exception("Not implemented (C)")
        # Cl included:
        if self.n == 6 and self.circular:
            self.pointgroup = C2v(n=self.n + len(self.bound_cl_to_c_positions))
        else:
            raise Exception("Not implemented (C-Cl)")

    def set_to_s_orbitals(self):
        # can't be undone!!!
        self.s_orbital_active = True
        for s in range(len(self.pointgroup.operations)):
            old_transform = self.pointgroup.operations[s].transform_p
            self.pointgroup.operations[s].transform_p = lambda x, f=old_transform: abs(f(x))


    def get_reducible_representation_for_ring_p_orbitals(self):
        reducible_representation = {}
        for op in self.pointgroup.operations:
            sum = 0
            for orbital in range(1, self.n+len(self.bound_cl_to_c_positions)+1):
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
        if self.s_orbital_active:
            p_symbols = sp.symbols(f"s1:{self.n +len(self.bound_cl_to_c_positions) + 1}")
        else:
            p_symbols = sp.symbols(f"p1:{self.n +len(self.bound_cl_to_c_positions) + 1}")
        SALCs = []
        for irred_name, value in irreducible_representation.items():
            # print("IRRED", irred_name)
            eq_expr = 0
            collected_p_orbitals_in_expression = {}
            if value != 0:
                irred = self.get_irreducible_representation_by_name(irred_name)
                for op in self.pointgroup.operations:
                    op_name = op.name
                    op_value = irred.characters[op_name] if op_name in irred.characters.keys() else 0
                    if op_value != 0:
                        if op.amount != 1:
                            raise Exception("Cannot handle this anymore")
                        transformed_p = op.transform_p(p_orbital_index)
                        # print("\t", p_orbital_index, "->", transformed_p)
                        coeff = op_value if transformed_p > 0 else -op_value
                        transformed_p = p_symbols[abs(transformed_p) - 1]
                        eq_expr += coeff * transformed_p  #* irred.dimension
                        collected_p_orbitals_in_expression[transformed_p] = collected_p_orbitals_in_expression.get(transformed_p, 0) + coeff
                        # print("\tadding", op.name, "=", coeff, f"{transformed_p}")

                eq_expr /= self.pointgroup.group_order()
                collected_p_orbitals_in_expression = {k: Fraction(v, self.pointgroup.group_order()) for k, v in collected_p_orbitals_in_expression.items()}
                if eq_expr != 0:
                    s = SALC(n=self.n+len(self.bound_cl_to_c_positions),
                         irred=irred_name, p_orbital_prefactors=collected_p_orbitals_in_expression,
                         equation=eq_expr, orbital_symbol="s" if self.s_orbital_active else "p")
                    # s.print()
                    SALCs.append(s)
        return SALCs


    def get_linear_independent_SALCs(self, salc_list:list, irreducible_representation):
        if linear_independent(salc_list):
            return salc_list
        # linear dependent -> orthogonalization:#TODO
        dict_of_salcs_grouped_by_irred = defaultdict(list)
        for salc in salc_list:
            dict_of_salcs_grouped_by_irred[salc.irred].append(salc)
        expected_numbers = self.get_expected_number_of_SALC_per_irreducible_representation(irreducible_representation)
        new_salc_list = []
        for key, items in dict_of_salcs_grouped_by_irred.items():
            if expected_numbers[key] == len(items):
                for i in items:
                    new_salc_list.append(i)
            elif expected_numbers[key] > len(items):
                raise Exception(f"Zu wenig SALCs für {key}")
            else:
                new_set = get_linear_independent_SALCs(items, expected_no=expected_numbers[key])
                # print(key, len(new_set))
                for s in new_set:
                    new_salc_list.append(s)
        return new_salc_list

    def get_all_SALCs(self):
        reducible_representation = self.get_reducible_representation_for_ring_p_orbitals()
        # print(reducible_representation)

        irreducible_representation = self.decomposing_into_irreducible_representations(reducible_representation)
        # print(irreducible_representation)

        all_SALCs = []
        for n in range(1, self.n +len(self.bound_cl_to_c_positions) + 1):
            SALCs = self.project(irreducible_representation, p_orbital_index=n)
            for s in SALCs:
                # s.print()
                checking_list = [is_multiple(x.equation, s.equation) for x in all_SALCs]# negative version of LC is still the same linear combination
                if not True in checking_list:
                    all_SALCs.append(s)

        # print(len(all_SALCs), irreducible_representation, sum(irreducible_representation.values()))
        all_SALCs = self.get_linear_independent_SALCs(all_SALCs, irreducible_representation)


        # expectation:
        # print( len(all_SALCs) , sum(self.get_expected_number_of_SALC_per_irreducible_representation(irreducible_representation).values()))
        assert len(all_SALCs) == sum(self.get_expected_number_of_SALC_per_irreducible_representation(irreducible_representation).values())
        return all_SALCs



    def orbitals_adjoint(self, p1, p2):
        alpha, beta = sp.symbols("alpha_s beta_s") if self.s_orbital_active else sp.symbols("alpha beta")
        if p1 == p2:
            return alpha
        if (p1 + 1 == p2 or p1 -1 == p2 ) and (p1 <= self.n and p2 <= self.n):
            return beta
        if self.circular and (
                            (p1-1 == 0 and p2 == self.n) or
                            (p2-1 == 0 and p1 == self.n)  ):
            return beta
        if len(self.bound_cl_to_c_positions) != 0 :
            if p1 in self.bound_cl_to_c_positions and p2 not in range(1,self.n+1):
                return beta
            if p2 in self.bound_cl_to_c_positions and p1 not in range(1,self.n+1):
                return beta
        return 0

    def h_eff(self, salc_1: SALC, salc_2: SALC):
        # print("h_eff:\t(", salc_1.equation, ") * (", salc_2.equation, end=")\t")
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
                        # if self.s_orbital_active:
                        #     print("s", p_orbital, "|", "s", combi, "=")
                        # else:
                        #     print("p", p_orbital, "|", "p", combi, "=")
                        # print("+", self.orbitals_adjoint(p_orbital, combi) ,"*", factor)
                h_eff_integral += part
        # print("=", h_eff_integral)
        return h_eff_integral

    def get_effective_hamilton_matrix(self, SALCs: list[SALC], irred:str):
        # print("get_effective_hamilton_matrix", [s.irred for s in SALCs], flush=True)
        # print(f"\n=== {irred} ===", type(irred), flush=True)

        n = len(SALCs)
        H = sp.zeros(n)
        H_show = sp.zeros(n)
        if len(SALCs) == 1:
                # h = sp.Symbol(f"H_{irred}")
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
        # for s in SALCs:
        #     s.print()
        SALCs_by_irred = norm_and_group_SALCs(SALCs)

        result = []
        alpha, beta, alpha_s, beta_s = sp.symbols(f"alpha beta alpha_s beta_s")
        sorting_dict_values = {alpha: 0, beta: -1, alpha_s: 0, beta_s: -1}
        for irred, salcs in SALCs_by_irred.items():
            H = self.get_effective_hamilton_matrix(SALCs=salcs, irred=irred)
            result_irred = calculate_hueckel_secular_equation(H, info=irred, sorting_dict_values=sorting_dict_values)
            for i in result_irred:
                i.symmetry = irred
                # print(i.symmetry)
                # sp.pprint(i.eigenvector)
                for row in range(len(i.eigenvector)):
                    if i.eigenvector[row] != 0:
                        i.salcs.append({"factor": i.eigenvector[row], "salc": salcs[row]})
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
