import sympy as sp

class SALC:
    def __init__(self, n:int, irred:str, p_orbital_prefactors:dict, equation:sp.Expr=None):
        self.n = n
        self.irred = irred
        self.equation = equation
        self.normalized = False

        self.prefactors_of_AOs = []
        for i in range(1, self.n+1):
            p_i = sp.symbols(f"p{i}")
            if p_i in p_orbital_prefactors.keys():
                self.prefactors_of_AOs.append( p_orbital_prefactors[p_i] )
            else:
                self.prefactors_of_AOs.append(0)
        # print("prefactors_of_AOs", self.prefactors_of_AOs)
        self.set_equation()

    def set_equation(self):
        if self.equation is not None:
            return
        p_symbols = sp.symbols(f"p1:{self.n + 1}")
        eq = 0
        for p in range(len(self.prefactors_of_AOs)):
            eq += p_symbols[p] * self.prefactors_of_AOs[p]
        self.equation = eq

    def print(self):
        print("SALC of", self.irred)
        sp.pprint(self.equation)
        print(self.prefactors_of_AOs)
        print()


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
            p_symbols = sp.symbols(f"p1:{self.n + 1}")
            self.equation = sum(coeff * p_symbols[idx]
                     for idx, coeff in enumerate(self.prefactors_of_AOs, 0))


        self.normalized = True
        self.print()


def norm_and_group_SALCs(list_of_SALCs:list[SALC]):
        SALCs_by_irred = {}
        for s in list_of_SALCs:
            s.norm()
            SALCs_by_irred.setdefault(s.irred, []).append(s)

        return SALCs_by_irred

