
import sympy as sp

from info_document.calculate_with_timeout import calculate_with_timeout
from molecule_orbital import molecule_orbital
from solve_2x2 import solve_2x2


def alternative_solver(H, info:str):
    S = sp.eye(H.rows)
    E = sp.symbols("E")
    secular_matrix = H - E * S

    det = sp.simplify(secular_matrix.det())
    # print("\n", info, ":")
    # sp.pprint(secular_matrix)
    # print("secular determinant:")
    # sp.pprint(det)

    # Solve the characteristic polynomial
    eigenvalues = sp.solve(det, E)

    # print("Eigenvalues:", len(eigenvalues))

    results = []
    for eigval in eigenvalues:
        results.append((eigval, None, []))

    return results


def calculate_hueckel_secular_equation(H, info: str, sorting_dict_values: dict):
    if H.rows != H.cols:
        raise Exception("Hamilton matrix has wrong dimensions")

    # Eigenwerte und Eigenvektoren bestimmen
    # H=H.subs(sp.Symbol('alpha_Cl'), sp.Symbol('alpha')*sp.Symbol('x'))
    # H=H.subs(sp.Symbol('beta_Cl'), sp.Symbol('beta') * sp.Symbol('x'))
    eigen_data = calculate_with_timeout(H.eigenvects, (), 60)
    # print("Eigen data:", "\n\t", eigen_data)
    if eigen_data is None:
        if H.rows == H.cols and H.rows == 2:
            eigen_data = solve_2x2(H=H)
        if eigen_data is None:
            eigen_data = alternative_solver(H=H, info=info)
    eigen_pairs = []

    for val, mult, vecs in eigen_data:
        for vec in vecs:
            eigen_pairs.append((val, vec.normalized()))

    # Sortieren nach substituierten Werten
    try:
        # print(eigen_pairs)
        # print(sorting_dict_values)
        sorted_pairs = sorted(
            eigen_pairs,
            key=lambda x: float( x[0].subs(sorting_dict_values) )
        )
    except TypeError:
        # cannot determine truth value of Relational:
        print("\tunsorted eigen pairs")
        sorted_pairs = eigen_pairs

    # print("\nEigenvalues (including degenerate orbitals):")
    molecule_orbitals = []
    for val, vec in sorted_pairs:
        m = molecule_orbital()
        m.eigenvalue = val
        if val == 0:
            raise Exception("0 eigenvalue")
        m.eigenvector = vec
        molecule_orbitals.append(m)
        # print("E =", end=" ")
        # sp.pprint(val)
        # print("v =", end=" ")
        # sp.pprint(vec)
        # print()
    return molecule_orbitals






