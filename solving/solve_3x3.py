
import sympy as sp

from solve_for_eigenvalues_first import find_eigenvector_general


def solve_3x3(H:sp.Matrix):
    if H.rows != 3 or H.cols != 3:
        raise Exception("this method is for 3x3 matrices only")
    # Eigen-Values:
    a = H[0, 0]
    b = H[0, 1]
    c = H[0, 2]
    d = H[1, 0]
    e = H[1, 1]
    f = H[1, 2]
    g = H[2, 0]
    h = H[2, 1]
    i = H[2, 2]
    assert b == d
    assert c == g
    assert f == h

    # j = -e-a-i
    # k = +a*i + a*e+e*i - f*f - b*b-c*c
    # l = -a*e*i + a*f*f + b*b*i - 2*c*b*f+ c*c*e
    j,k,l = sp.symbols("j k l")
    x = sp.symbols("x")
    eq = x ** 3 + j * x ** 2 + k * x + l
    roots = sp.solve(eq, x)

    mults = []
    epsilons = []
    for root in roots:
        mults.append(None)
        epsilons.append(root)

    # Eigen-Vectors:
    eigenvecs = []
    for eigenvalue in epsilons:
        vector = find_eigenvector_general(H, eigenvalue=eigenvalue)
        eigenvecs.append(vector)

    # collect:
    if len(mults) != len(epsilons) or len(mults) != len(eigenvecs):
        raise Exception("strange")
    eigen_data = []
    for i in range(len(mults)):
        eigen_data.append([epsilons[i], mults[i], eigenvecs[i]])
    return eigen_data



# alpha, beta, beta_Cl, alpha_Cl = sp.symbols("alpha beta beta_Cl alpha_Cl")
# H = sp.Matrix([
#     [alpha, sp.sqrt(2)*beta, 2*beta_Cl],
#     [sp.sqrt(2)*beta, alpha + beta, 0],
#     [2*beta_Cl, 0, alpha_Cl]
# ])

# Temporäre Symbole
a, b, c, d, e, f, g, h, i = sp.symbols("a b c d e f g h i")
H = sp.Matrix([
    [a, b, c],
    [b, e, f],
    [c, f, i]
])

x = solve_2x2(H)
print(x)