


import sympy as sp

def solve_for_eigenvalues_first(H: sp.Matrix):
    """
    Returns eigenvalue data in the same structure as H.eigenvects(),
    but without computing eigenvectors.
    """
    print("solve_for_eigenvalues_only", H, "...", flush=True)
    epsilons = []
    mults = []
    for value, mult in H.eigenvals().items():
        epsilons.append(value)
        mults.append(mult)

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


def find_eigenvector_general(H: sp.Matrix, eigenvalue):
    """
        Returns a list of eigenvectors corresponding to the given eigenvalue.
        Each eigenvector is a sympy Matrix.
    """
    M = H - eigenvalue * sp.eye(H.shape[0])
    print("find_eigenvector_general", M, "\n", eigenvalue, "...", flush=True)
    nullspace = M.nullspace()
    return nullspace # list of vectors


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