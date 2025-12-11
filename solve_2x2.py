

import sympy as sp



def solve_2x2(H:sp.Matrix):
    if H.rows != 2 or H.cols != 2:
        raise Exception("this method is for 2x2 matrices only")
    # Eigen-Values:
    a = H[0, 0]
    b = H[0, 1]
    c = H[1, 0]
    d = H[1, 1]

    p = -(a + d)
    q = a * d - c * b
    epsilon_1 = - (p / 2) + sp.sqrt((p / 2) ** 2 - q)
    epsilon_2 = - (p / 2) - sp.sqrt((p / 2) ** 2 - q)

    # Eigen-Vectors:
    vec_1 = find_eigenvector_2x2(H, eigenvalue=epsilon_1)
    vec_2 = find_eigenvector_2x2(H, eigenvalue=epsilon_2)

    return [epsilon_1, 1, vec_1], [epsilon_2, 1, vec_2]


def find_eigenvector_2x2(H: sp.Matrix, eigenvalue):
    x, y = sp.symbols("x y")
    a = H[0, 0]
    b = H[0, 1]
    c = H[1, 0]
    d = H[1, 1]

    eqs = [
        (a - eigenvalue)*x + b*y,
        c*x + (d - eigenvalue)*y
    ]

    sol = sp.solve(eqs, (x, y), dict=True)

    if len(sol) != 1:
        # Regelfall
        sol = sp.solve([eqs[0].subs(x, 1), eqs[1].subs(x, 1)], (y,), dict=True)
        return sp.Matrix([1, sol[0][y]])

    print("strange!")
    return sp.Matrix([sol[0][x], sol[0][y]])



if __name__ == '__main__':
    a,b,c,d, p, q = sp.symbols("a b c d p q")
    H = sp.Matrix([[a, b],[c,d]])

    eigen_data = solve_2x2(H)

    for val, mult, vecs in eigen_data:
        sp.pprint(vecs)