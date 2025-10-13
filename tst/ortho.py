import sympy as sp

from SALC import SALC, get_linear_independent_SALCs





if __name__ == "__main__":
    irred = "E1g"
    p1, p2, p3, p4, p5, p6 = sp.symbols(f"p1 p2 p3 p4 p5 p6")
    s1 = SALC(n=6, p_orbital_prefactors={p1: 1/6, p2: +2/12, p3: -1/12, p4: -1/6, p5: -1/12, p6: 1/12},
              irred=irred, orbital_symbol="p")
    s2 = SALC(n=6, p_orbital_prefactors={p1: 1 / 12, p2: +1/ 6, p3: 1/ 12, p4: -1 / 12, p5: -1 / 6, p6: -1 / 12},
              irred=irred, orbital_symbol="p")
    s3 = SALC(n=6, p_orbital_prefactors={p1: -1 / 12, p2: +1 / 12, p3: 1 / 6, p4: 1 / 12, p5: -1 / 12, p6: -1 / 6},
              irred=irred, orbital_symbol="p")

    # print([s3.prefactors_of_AOs[i] + s2.prefactors_of_AOs[i] for i in range(6)])
    x = get_linear_independent_SALCs([s1, s2, s3], expected_no = 2)
    for s in x:
        s.norm()
        s.print()
    # print(len(x))
    assert len(x) == 2

