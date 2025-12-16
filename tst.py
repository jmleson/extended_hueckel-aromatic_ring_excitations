import sympy as sp

from solve_for_eigenvalues_first import solve_for_eigenvalues_first, find_eigenvector_general


def solve_with_temp_symbols(H, subs_in):
    """
    Solve eigenvalues using solve_for_eigenvalues_only by replacing complex expressions
    with temporary symbols.
    """
    # 1. Ersetze komplizierte Ausdrücke durch temporäre Symbole
    H_subs = H.subs(subs_in)

    # 2. Eigenwerte berechnen
    eigen_data_temp = solve_for_eigenvalues_first(H_subs)

    # 3. Original-Ausdrücke wieder einsetzen
    subs_out = {v: k for k, v in subs_in.items()}
    eigen_data_original = []
    for val, mult, vectors in eigen_data_temp:
        val_orig = val.subs(subs_out)
        vectors_orig = [v.subs(subs_out) for v in vectors]
        eigen_data_original.append([val_orig, mult, vectors_orig])

    return eigen_data_original


# -----------------------------
# Beispiel für deine Matrix
alpha, beta, beta_Cl, alpha_Cl = sp.symbols("alpha beta beta_Cl alpha_Cl")

H = sp.Matrix([
    [alpha, sp.sqrt(2)*beta, 2*beta_Cl],
    [sp.sqrt(2)*beta, alpha + beta, 0],
    [2*beta_Cl, 0, alpha_Cl]
])

# Temporäre Symbole
a, b, c, d, e, f, g, h, i = sp.symbols("a b c d e f g h i")
H = sp.Matrix([
    [a, b, c],
    [b, e, f],
    [c, f, i]
])

x = solve_for_eigenvalues_first(H)
print(x)

x = find_eigenvector_general(H, eigenvalue = sp.Symbol("g"))
print(x)


# x = sp.symbols("x")
# charpoly = H.charpoly(x)
# roots = sp.solve(charpoly.as_expr(), x)  # kann CRootOf enthalten
# for r in roots:
#     print(r)