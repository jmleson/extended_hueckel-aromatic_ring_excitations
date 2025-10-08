import sympy as sp

from calculate_hueckel_secular_equation import calculate_hueckel_secular_equation
from molecule_orbital import molecule_orbital

#
# def calculate_hueckel_secular_equation(H, info: str, sorting_dict_values: dict):
#     if H.rows != H.cols:
#         raise Exception("Hamilton matrix has wrong dimensions")
#
#     S = sp.eye(H.rows)
#     E = sp.symbols("E")
#     secular_matrix = H - E * S
#
#     # det = sp.simplify(secular_matrix.det())
#     # print("\n", info, ":")
#     # print("secular determinant:")
#     # sp.pprint(det)
#
#     # Eigenwerte und Eigenvektoren bestimmen
#     eigen_data = H.eigenvects()
#     eigen_pairs = []
#
#     for val, mult, vecs in eigen_data:
#         for vec in vecs:
#             eigen_pairs.append((val, vec.normalized()))
#
#     # Sortieren nach substituierten Werten
#     sorted_pairs = sorted(
#         eigen_pairs,
#         key=lambda x: x[0].subs(sorting_dict_values)
#     )
#
#     # print("\nEigenvalues (including degenerate orbitals):")
#     molecule_orbitals = []
#     for val, vec in sorted_pairs:
#         m = molecule_orbital()
#         m.eigenvalue = val
#         m.eigenvector = vec
#         molecule_orbitals.append(m)
#         # print("E =", end=" ")
#         # sp.pprint(val)
#         # print("v =", end=" ")
#         # sp.pprint(vec)
#         # print()
#     return molecule_orbitals


def get_six_orbitals_in_benzene_form(alpha_symbol, beta_symbol, info:str="6 orbitals"):
    alpha, beta = sp.symbols(f'{alpha_symbol} {beta_symbol}')
    H = sp.Matrix([
        [alpha, beta,  0,     0,     0,     beta ],
        [beta,  alpha, beta,  0,     0,     0    ],
        [0,     beta,  alpha, beta,  0,     0    ],
        [0,     0,     beta,  alpha, beta,  0    ],
        [0,     0,     0,     beta,  alpha, beta ],
        [beta,  0,     0,     0,     beta,  alpha ]
    ])
    return calculate_hueckel_secular_equation(H, info=info, sorting_dict_values = {alpha: 0, beta: -1})


def get_benzene():
    p_orbitals = get_six_orbitals_in_benzene_form(alpha_symbol="alpha", beta_symbol="beta", info="p orbitals")
    s_orbitals = get_six_orbitals_in_benzene_form(alpha_symbol="alpha_s", beta_symbol="beta_s", info="s orbitals")
    return {"p": p_orbitals, "s": s_orbitals}


if __name__ == "__main__":
    get_benzene()





