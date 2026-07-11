


import sympy as sp

from info_document.calculate_with_timeout import calculate_with_timeout
from solving.old.solve_2x2 import find_eigenvector_2x2


def solve_for_eigenvalues_first(H: sp.Matrix):
    """
    Returns eigenvalue data in the same structure as H.eigenvects(),
    but without computing eigenvectors.
    """
    print("solve_for_eigenvalues_only", "...", flush=True)
    epsilons = []
    mults = []
    try:
        for value, mult in H.eigenvals().items():
            epsilons.append(value)
            mults.append(mult)
    except Exception as exc:
        print(exc, flush=True)
        return []

    # Eigen-Vectors:
    eigenvecs = []
    for eigenvalue in epsilons:
        vector = find_eigenvector_general(H, eigenvalue=eigenvalue)
        # print("type(vector)", type(vector), flush=True)
        eigenvecs.append(vector)

    # collect:
    if len(mults) != len(epsilons) or len(mults) != len(eigenvecs):
        raise Exception("strange")
    eigen_data = []
    for i in range(len(mults)):
        eigen_data.append([epsilons[i], mults[i], eigenvecs[i]])
    return eigen_data


def find_eigenvector_general(H: sp.Matrix, eigenvalue) -> list[sp.Matrix]:
    """
        Returns a list of eigenvectors corresponding to the given eigenvalue.
        Each eigenvector is a sympy Matrix.
    """
    M = H - eigenvalue * sp.eye(H.shape[0])
    # print("find_eigenvector_general", "...", flush=True)
    nullspace = calculate_with_timeout(M.nullspace, (), timeout_in_s=5*60, msg=f"Calculation timed out after {10} s for M.nullspace()")
    # print("\talready found?", nullspace is None, flush=True)
    if nullspace is None:
        nullspace = find_eigenvector_2x2(H=H, eigenvalue=eigenvalue)
        # print("\tafter find2x2", nullspace is None, flush=True)
    if nullspace is None:
        nullspace = find_3x3_eigenvectors(H, eigenvalue=eigenvalue)
        # print("\tafter find3x3", nullspace is None, flush=True)
    if nullspace is None:
        nullspace = construct_empty_eigenvector(H, eigenvalue)
        # print("\tafter construct_empty_eigenvector", nullspace is None, flush=True)
    return nullspace # list of vectors


def find_3x3_eigenvectors(H: sp.Matrix, eigenvalue):
    if H.rows != 3 or H.cols != 3:
        # raise Exception("this method is for 3x3 matrices only")
        return None
    M = H - eigenvalue * sp.eye(H.shape[0])
    if not M.is_symmetric():
        return None
    a = H[0, 0]
    b = H[0, 1]
    c = H[0, 2]
    e = H[1, 1]
    f = H[1, 2]

    v1 = (-b*f+c*(e - eigenvalue))/(b**2 - (a-eigenvalue)*(e-eigenvalue))
    v2 = (-a*f+b*c+eigenvalue*f)/(a*e-a*eigenvalue-b*b -e*eigenvalue+eigenvalue**2)
    v3 = 1
    vec = sp.Matrix([[v1], [v2], [v3]])
    return [vec]


def construct_empty_eigenvector(H:sp.Matrix, eigenvalue):
    return [sp.Matrix([sp.symbols(f'k{i}') for i in range(H.rows)])]



if __name__ == '__main__':
    # Symbole definieren
    alpha, beta = sp.symbols("alpha beta")
    beta_Cl, alpha_Cl = sp.symbols("beta_Cl alpha_Cl")

    # Chlorbenzen-Problemmatrix:
    H = sp.Matrix([
        [alpha, sp.sqrt(2) * beta, 0, 0, beta_Cl],
        [sp.sqrt(2) * beta, alpha, beta, 0, 0],
        [0, beta, alpha, sp.sqrt(2) * beta, 0],
        [0, 0, sp.sqrt(2) * beta, alpha, 0],
        [beta_Cl, 0, 0, 0, alpha_Cl]
    ])
    # C6H4Cl2-Problemmatrix:
    H = sp.Matrix([[alpha, sp.sqrt(2)*beta, 2*beta_Cl],
                   [sp.sqrt(2)*beta, alpha + beta, 0],
                   [2*beta_Cl, 0, alpha_Cl]])

    x = solve_for_eigenvalues_first(H)
    print(x)