# import sympy as sp
#
# alpha, alpha_s, beta, beta_s, delta, A, B = sp.symbols('alpha alpha_s beta beta_s delta A B')
#
# E_SQ = ( 32 * delta**4 ) / (-2 * alpha + 2 * alpha_s - 2 * beta + 2 * beta_s) + (32 * delta**4)/(-2*alpha + 2*alpha_s)
# # sp.pprint(E_SQ)
#
# E_TT = (4 * delta**4)/(2 * alpha + 2 * alpha_s + 2 * beta - 2 * beta_s) + (36 * delta**4)/(-2*alpha + 2*alpha_s -2*beta+2*beta_s) + (24*delta**4)/(-2*alpha + 2*alpha_s)
#
# E_SS = (32 * delta**4)/(-alpha + alpha_s - beta + beta_s)
#
# eq1 = sp.Eq(A, alpha - alpha_s)
# eq2 = sp.Eq(B, beta - beta_s)
#
#
#
# solution = sp.solve([E_SQ, E_TT, E_SS, eq1, eq2], (A, B))
# print(solution)

#
import sympy as sp

A, B, d = sp.symbols('A B delta')
E_SQ ,E_TT ,E_SS = sp.symbols('E_SQ E_TT E_SS')

eq1 = (32*d**4)/(2*A) + (32 *d**4)/(2*B)  -E_SQ
eq2 = (4*d**4)/(4*B-2*A) + (36*d**4)/(2*A) + (24*d**4)/(2*B)   - E_TT
eq3 = (32*d**4)/A  - E_SS

# Gleichungssystem nach A und B lösen
solution = sp.solve([eq1, eq2, eq3], (A, B, d), dict=True)
print("Symbolische Lösungen:")
print(solution)
