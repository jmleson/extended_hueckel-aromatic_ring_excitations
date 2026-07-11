
import sympy as sp

def get_eigenvalues(H: sp.Matrix, info:str) -> tuple[list[sp.Expr], list[int]]:
    """
    find eigenvalues (energies of MOs)
    """
    # if H.rows != 2 or H.cols != 2:
    #     print("TODO 2x2 solver?? ")

    epsilons = []
    mults = []
    symbols = []
    try:
        H.eigenvals()
        for value, mult in H.eigenvals().items():
            epsilons.append(value)
            mults.append(mult)
    except Exception as exc:
        print(exc)
        # un-solvable -> use general variables
        epsilons = [sp.Symbol(f"E{nr}({info})") for nr in range(H.rows)]
        symbols.extend(epsilons)
        mults = [1 for nr in range(H.rows)]

    assert len(mults) == len(epsilons)
    return epsilons, mults, symbols

