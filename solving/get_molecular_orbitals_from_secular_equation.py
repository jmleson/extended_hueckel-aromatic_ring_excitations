import sympy as sp

from molecule_orbital import molecule_orbital
from solving.get_eigenvalues import get_eigenvalues
from solving.get_eigenvectors_from_eigenvalue import get_eigenvectors_from_eigenvalue


def get_molecular_orbitals_from_secular_equation(H: sp.Matrix, info:str, sorting_dict_values:dict):
    if H.rows != H.cols:
        raise Exception("Hamilton matrix has wrong dimensions")

    epsilons, mults, symbols = get_eigenvalues(H=H, info=info)
    for symbol in symbols:
        sorting_dict_values[symbol] = 10

    eigen_vectors = []
    eigen_values = []
    for i in range(len(epsilons)):
        _eigen_vectors_ = get_eigenvectors_from_eigenvalue(H=H, eigenvalue=epsilons[i], expected_mult=mults[i])
        assert len(_eigen_vectors_) == mults[i]
        for v in _eigen_vectors_:
            eigen_vectors.append(v)
            eigen_values.append(epsilons[i])

    # COLLECT DATA INTO MOLECULAR_ORBITALS:
    # if not sum(mults) == len(eigen_vectors):
    #     print(mults)
    #     print(eigen_vectors)
    assert sum(mults) == len(eigen_vectors)

    molecule_orbitals = []
    for i in range(sum(mults)):
        m = molecule_orbital()
        m.eigenvalue = eigen_values[i]
        # if m.eigenvalue == 0:# TODO reactivate
        #     raise Exception("0 eigenvalue")
        m.eigenvector = eigen_vectors[i]
        molecule_orbitals.append(m)
        # print("E =", end=" ")
        # sp.pprint(val)
        # print("v =", end=" ")
        # sp.pprint(vec)
        # print()

    # SORT ORBITALS ACCORDING TO ENERGY:
    try:
        sorted_pairs = sorted(
                    molecule_orbitals,
                    key=lambda x: float( x.eigenvalue.subs(sorting_dict_values) )
                )
        print("\t", len(sorted_pairs), "MOs for ", H.rows, "x", H.cols, "matrix, mult-sum = ", sum(mults))
    except:
        sorted_pairs = molecule_orbitals
    assert len(sorted_pairs) == H.rows or len(sorted_pairs) == sum(mults)
    return sorted_pairs




