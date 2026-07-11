import sympy as sp

from is_multiple import is_multiple


class SALC:
    def __init__(self, n:int, irred:str, p_orbital_prefactors:dict, orbital_symbol:str, equation:sp.Expr=None):
        self.n = n
        self.irred = irred
        self.equation = equation
        self.normalized = False
        self.orbital_symbol = orbital_symbol

        self.prefactors_of_AOs = []
        for i in range(1, self.n+1):
            p_i = self.get_symbols()[i-1]
            if p_i in p_orbital_prefactors.keys():
                self.prefactors_of_AOs.append( p_orbital_prefactors[p_i] )
            else:
                self.prefactors_of_AOs.append(0)
        # print("prefactors_of_AOs", self.prefactors_of_AOs)
        self.set_equation()

    def get_symbols(self):
        return sp.symbols(f"{self.orbital_symbol}1:{self.n + 1}")

    def set_equation(self):
        if self.equation is not None:
            return
        p_symbols = self.get_symbols()
        eq = 0
        for p in range(len(self.prefactors_of_AOs)):
            eq += p_symbols[p] * self.prefactors_of_AOs[p]
        self.equation = eq

    def print(self):
        print("SALC of", self.irred)
        sp.pprint(self.equation)
        print(self.prefactors_of_AOs)
        print()

    def multiple_in_list(self, list_of_other_salcs:list):
        for i in list_of_other_salcs:
            if is_multiple(i.equation, self.equation):
                return True
        return False


    def norm(self):
        if self.normalized:
            return
        # assumes, that included p orbitals are normalized
        # self.print()
        current_norm_value = sum(coeff ** 2 for coeff in self.prefactors_of_AOs)
        # print(current_norm_value)

        if current_norm_value != 1:
            self.prefactors_of_AOs = [x/sp.sqrt(current_norm_value) for x in self.prefactors_of_AOs]

            #update equation:
            p_symbols = sp.symbols(f"{self.orbital_symbol}1:{self.n + 1}")
            self.equation = sum(coeff * p_symbols[idx]
                     for idx, coeff in enumerate(self.prefactors_of_AOs, 0))


        self.normalized = True
        # self.print()


def norm_and_group_SALCs(list_of_SALCs:list[SALC]):
        SALCs_by_irred = {}
        for s in list_of_SALCs:
            s.norm()
            if s.irred not in SALCs_by_irred.keys():
                SALCs_by_irred[s.irred] = [s]
            else:
                SALCs_by_irred[s.irred] += [s]
        return SALCs_by_irred


def linear_independent( salc_list:list):
    A = sp.Matrix([
        [sp.Rational(val) if not isinstance(val, sp.Basic) else val
         for val in salc.prefactors_of_AOs]
        for salc in salc_list
    ])
    rank = A.rank()
    #rank == len(salc_list):
    # not linearly dependent
    return rank == len(salc_list)


def get_linear_independent_SALCs(salc_list: list[SALC], expected_no:int) -> list[SALC]:
    if len(salc_list) == expected_no:
        return salc_list
    if linear_independent(salc_list):
        return salc_list
    if not salc_list:
        return []
    set_of_irreds = set([i.irred for i in salc_list])
    if len(set_of_irreds) != 1:
        raise Exception("method supposed to use with one irred only")

    n = salc_list[0].n
    orbital_symbol = salc_list[0].orbital_symbol
    p = sp.symbols(f"{orbital_symbol}1:{n+1}")

    # take last two salcs and build linear combination:
    s1, s2 = salc_list[-2], salc_list[-1]
    new_prefactors = [s1.prefactors_of_AOs[i] + s2.prefactors_of_AOs[i] for i in range(n)]
    new_salc = SALC(
        n=n,
        p_orbital_prefactors={p[i]: new_prefactors[i] for i in range(n)},
        irred=s1.irred,
        orbital_symbol=orbital_symbol
    )
    if new_salc.multiple_in_list(salc_list):
        new_prefactors = [s1.prefactors_of_AOs[i] - s2.prefactors_of_AOs[i] for i in range(n)]
        new_salc = SALC(
            n=n,
            p_orbital_prefactors={p[i]: new_prefactors[i] for i in range(n)},
            irred=s1.irred,
            orbital_symbol=orbital_symbol
        )
        # 2nd test:
        if new_salc.multiple_in_list(salc_list):
            raise Exception("case not implemented yet")
    new_list = salc_list[:-2] + [new_salc]
    return get_linear_independent_SALCs(new_list, expected_no)



# def combine_eigenvectors_to_linearily_combined_salcs_analogously(eigenvectors:list, expected_no):
#     # TODO



if __name__ == "__main__":
    p1, p2, p3, p4, = sp.symbols(f"p1 p2 p3 p4")
    s = SALC(n=4, irred="A2u",
         p_orbital_prefactors={p1: 1 / sp.sqrt(4), p2: 1 / sp.sqrt(4), p3: 1 / sp.sqrt(4),
                               p4: 1 / sp.sqrt(4)},
         orbital_symbol="p")

    print(s.prefactors_of_AOs)