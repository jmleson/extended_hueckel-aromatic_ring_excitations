# SALCs for benzene (n=6) with exact symbolic arithmetic (Rationals and roots)
import sympy as sp

def benzene_salcs_symbolic(n=6):
    # assert n == 6, "This routine is written for the 6-membered ring (benzene)."
    # Basis symbols phi1..phi6
    phi = sp.symbols(f'phi1:{n+1}')  # phi1 ... phi6

    # root of unity omega = exp(2*pi*i/n)
    k_sym = sp.symbols('k', integer=True)
    omega = sp.exp(2*sp.pi*sp.I / n)

    # build complex eigenvectors v_k with components omega^(k*j), j=0..n-1
    def v_k(k):
        return sp.Matrix([omega**(k * j) for j in range(n)])

    # generate v_k for k=0..n-1
    vs = [v_k(k) for k in range(n)]

    # Convert complex vectors into real SALCs:
    # - for k=0 and k=n/2 (here 3) use the vector itself (real)
    # - for 1<=k<=n/2-1 produce two real vectors: Re(v_k) and Im(v_k)
    real_vecs = []

    # k = 0
    real_vecs.append(sp.simplify(sp.Matrix([sp.re(c) for c in vs[0]])))

    # k = 1 and k = 5 => pair -> Re and Im
    for k in (1, 2):
        vk = vs[k]
        v_minus = vs[n-k]
        # real part: (v_k + v_{-k})/2  -> equals Re(v_k)
        re_vk = sp.simplify((vk + v_minus) / 2)
        # imag part: (v_k - v_{-k})/(2I) -> equals Im(v_k)
        im_vk = sp.simplify((vk - v_minus) / (2*sp.I))
        # convert components to exact trig forms (cos/sin where possible)
        # but SymPy normally keeps trig exact for rational multiples of pi
        re_vk = sp.Matrix([sp.simplify(sp.trigsimp(sp.re(c))) for c in re_vk])
        im_vk = sp.Matrix([sp.simplify(sp.trigsimp(sp.re(c))) for c in im_vk])  # re() to remove tiny imag=0
        real_vecs.append(re_vk)
        real_vecs.append(im_vk)

    # k = 3 (n/2)
    real_vecs.append(sp.simplify(sp.Matrix([sp.re(c) for c in vs[3]])))

    # Now real_vecs is a list of 6 real vectors (symbolic)
    # Express each SALC as combination of phi_j
    salcs = []
    for vec in real_vecs:
        # components can contain cos(...) etc. we simplify them
        comps = [sp.simplify(sp.trigsimp(sp.N(c, n=50) if c.has(sp.I) else c)) for c in vec]
        # better to convert trig to exact rationals/sqrts if possible:
        comps = [sp.simplify(sp.nsimplify(c, [sp.sqrt(3)])) for c in vec]
        expr = sum(comps[j] * phi[j] for j in range(n))
        expr = sp.simplify(expr)
        salcs.append(expr)

    # Normalize each SALC symbolically: compute norm = sqrt(sum(c^2))
    salcs_normed = []
    for expr in salcs:
        # extract coefficients for phi basis
        coeffs = [sp.simplify(sp.expand(expr).coeff(p)) for p in phi]
        # compute norm^2 symbolically
        norm2 = sp.simplify(sum(sp.simplify(sp.expand(c**2)) for c in coeffs))
        # norm might produce rational or expression with sqrt(3)
        norm = sp.sqrt(sp.simplify(norm2))
        # form normalized SALC (if norm==0 skip)
        if norm == 0:
            salcs_normed.append(expr)  # degenerate zero (shouldn't happen)
        else:
            salcs_normed.append(sp.simplify(expr / norm))

    return phi, salcs_normed

if __name__ == "__main__":
    phi, salcs = benzene_salcs_symbolic(4)
    print("Basis:", phi)
    print("\nSix SALCs (symbolically normalized):\n")
    for i, s in enumerate(salcs, start=1):
        print(f"SALC {i} =")
        sp.pprint(s)
        print()
