import sympy as sp

# Symbole
alpha, beta, beta_Cl, alpha_Cl = sp.symbols("alpha beta beta_Cl alpha_Cl")
x = sp.symbols("x")

a, b, c, d, e, f, g, h, i = sp.symbols("a b c d e f g h i")
H = sp.Matrix([
    [a, b, c],
    [b, e, f],
    [c, f, i]
])

# 1. Eigenwerte symbolisch
roots = sp.solve(H.charpoly(x).as_expr(), x)
lambda1, lambda2, lambda3 = roots

# 2. Eigenvektoren berechnen
eigenvectors = []
for lam in roots:
    v = (H - lam*sp.eye(3)).nullspace()
    eigenvectors.append(v)

for lam, v in zip(roots, eigenvectors):
    print("Eigenvalue:", lam)
    print("Eigenvector(s):", v)
