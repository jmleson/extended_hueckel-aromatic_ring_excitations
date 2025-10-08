import sympy as sp

def round_and_collect(expr, symbols, digits=3):
    """Rundet alle Zahlen in einem SymPy-Ausdruck und fasst nach den gegebenen Symbolen zusammen."""
    expr_approx = expr.evalf(digits + 2)  # kleine Reserve für Rundungsgenauigkeit
    rounded_expr = expr_approx.xreplace({
        n: round(float(n), digits)
        for n in expr_approx.atoms(sp.Float)
    })
    return sp.collect(sp.simplify(rounded_expr), symbols)
