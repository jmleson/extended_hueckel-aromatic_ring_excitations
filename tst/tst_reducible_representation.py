import sympy as sp

n = 6
orbitals = list(range(n))

# Orbitaltypen (z.B. 'pz' für π-Orbitale senkrecht zur Ebene)
orbital_type = ['pz']*n

# Symmetrieoperationen als Abbildung der Indizes und Vorzeichen
def E(i): return i, 1
def C6(i): return (i + 1) % n, 1
def C6_5(i): return (i + 5) % n, 1
def C3(i): return (i + 2) % n, 1
def C3_2(i): return (i + 4) % n, 1
def C2(i): return (i + 3) % n, 1
def C2_prime(i): return [0,3,4,1,2,5][i], 1
def C2_double_prime(i): return [0,3,5,2,4,1][i], 1

def sigma_h(i):
    # σh invertiert pz-Orbitale
    return (i + 3) % n, -1

def sigma_v(i): return [0,5,4,3,2,1][i], 1
def sigma_d(i): return [0,3,2,5,4,1][i], 1
def inversion(i):
    # pz-Orbitale ungerade unter i → -1
    return (i + 3) % n, -1

def S6(i):
    j, s = C6(i)
    _, s_h = sigma_h(j)
    return j, s * s_h

def S6_5(i):
    j, s = C6_5(i)
    _, s_h = sigma_h(j)
    return j, s * s_h

operations = [E, C6, C6_5, C3, C3_2, C2, C2_prime, C2_double_prime,
              sigma_h, sigma_v, sigma_d, inversion, S6, S6_5]
labels = ['E', 'C6', 'C6^5', 'C3', 'C3^2', 'C2', "C2'", "C2''",
          'σh', 'σv', 'σd', 'i', 'S6', 'S6^5']

# Berechnung der reduzierten Darstellung (Charakter) mit Vorzeichen
characters = []
for R in operations:
    trace = sum(s for i in orbitals if R(i)[0] == i or R(i)[0] == -i for _, s in [R(i)])
    characters.append(trace)

for label, char in zip(labels, characters):
    print(f"{label:>5}: χ_red = {char}")
