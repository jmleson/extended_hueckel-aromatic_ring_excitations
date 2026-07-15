import sympy as sp


def is_multiple(eq_1: sp.Expr, eq_2: sp.Expr):
    # calculate ratio:
    ratio = sp.simplify(eq_1 / eq_2)
    # determine whether ratio is just a numerical factor:
    return ratio.is_number