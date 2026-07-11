
import sympy as sp

from info_document.calculate_with_timeout import calculate_with_timeout


def get_eigenvectors_from_eigenvalue(H: sp.Matrix, eigenvalue, expected_mult:int):
    """
           Returns a list of eigenvectors corresponding to the given eigenvalue.
           Each eigenvector is a sympy Matrix.
    """
    M = H - eigenvalue * sp.eye(H.shape[0])
    # if H.rows != 2 or H.cols != 2:
    #     print("TODO 2x2 solver?? ")

    to_be_normed = True
    try:
        vectors = calculate_with_timeout(M.nullspace, (), timeout_in_s=5*60, msg=f"Calculation timed out after {10} s for M.nullspace()")
    except Exception as exc:
        print(exc)
        vectors = []
    if vectors is None:
        vectors = []

    if len(vectors) < expected_mult:
        for i in range(expected_mult):
            vectors.append( construct_empty_eigenvector(H=H) )
        to_be_normed = False
    if len(vectors) > expected_mult:
        raise Exception("too many eigenvectors found!")

    # TRY NORMALIZATION
    try:
        assert to_be_normed # use true norm only in cases where we found an eigenvector
        vectors = [v.normalized() for v in vectors]
    except Exception:
        factor = sp.Symbol("Norm")
        vectors = [factor * v for v in vectors]
    assert len(vectors) == expected_mult
    return vectors


def construct_empty_eigenvector(H:sp.Matrix):
    norm = sp.Symbol("Norm")
    dummy_values = sp.symbols('?_1:%d' % (H.cols + 1))
    v = norm * sp.Matrix(dummy_values)
    return v


