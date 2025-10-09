import sympy as sp
from molecule_orbital import molecule_orbital


class TransitionIntegral:
    def __init__(self, bra: molecule_orbital, ket:molecule_orbital):
        # assumes dipole-z-operator
        self.bra = bra
        self.ket = ket

    def print(self):
         integral = "< " + str(self.bra.get_salc_equation()) + " | q_z | " + str(self.ket.get_salc_equation()) + " >"
         print(integral)
         return integral

    def get_multipliers_as_sum(self, linear_combination_mo):
        multipliers = []
        for s in linear_combination_mo.salcs:
            p_symbols = s["salc"].get_symbols()
            # print("eq = ", s["salc"].equation)
            for i in range(len(s["salc"].prefactors_of_AOs)):
                factor = s["salc"].prefactors_of_AOs[i]
                total_factor = s["factor"] * factor
                if total_factor != 0:
                    multipliers.append({"total_factor": total_factor, "p_symbols": p_symbols[i], "p_number": i+1})
            # print("->", multipliers)
        return multipliers

    def multiply_out(self):
        bra_multipliers = self.get_multipliers_as_sum(self.bra)
        ket_multipliers = self.get_multipliers_as_sum(self.ket)

        result = 0
        delta = sp.Symbol("delta")
        for bra_eq in bra_multipliers:
            for ket_eq in ket_multipliers:
                factor = bra_eq["total_factor"] * ket_eq["total_factor"]
                if factor != 0:
                    # eq = bra_eq["p_symbols"] * ket_eq["p_symbols"]
                    # print("...", eq)
                    #  |q_z | s > = p (shape of p orbital by multiplication with dipole operator)
                    if bra_eq["p_number"] != ket_eq["p_number"]:
                        # orbitals located at different atoms -> 0 in hückel approximation
                        factor = 0
                    result += factor * delta
        # print("result:", result)
        return result


