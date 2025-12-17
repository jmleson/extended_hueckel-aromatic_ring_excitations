# from calculate_hueckel_secular_equation import calculate_hueckel_secular_equation
# import sympy as sp
#
#
# heterosymbol = "Cl"
# alpha, beta, alpha_s, beta_s = sp.symbols(f"alpha beta alpha_s beta_s")
# alpha_Cl, beta_Cl, alpha_s_Cl, beta_s_Cl = sp.symbols(f"alpha_{heterosymbol} beta_{heterosymbol} alpha_s_{heterosymbol} beta_s_{heterosymbol}")
# alpha_N, beta_N, alpha_s_N, beta_s_N = sp.symbols("alpha_N beta_N alpha_s_N beta_s_N")
# sorting_dict_values = {alpha: 0, beta: -1,
#                                alpha_s: 0, beta_s: -1,
#                                alpha_Cl: 0, beta_Cl: -1,
#                                alpha_s_Cl: 0, beta_s_Cl: -1,
#                                alpha_N: 0, beta_N: -1,
#                                alpha_s_N: 0, beta_s_N: -1
#                                }
#
# irred="Ag"
# alpha, beta = sp.symbols("alpha beta")
# beta_Cl, alpha_Cl = sp.symbols("beta_Cl alpha_Cl")
#
# H = sp.Matrix([[alpha, sp.sqrt(2)*beta, 2*beta_Cl],
#                    [sp.sqrt(2)*beta, alpha + beta, 0],
#                    [2*beta_Cl, 0, alpha_Cl]])
# # H = sp.Matrix([[alpha+beta, 0],
# #                 [0 , alpha - beta],
# #                 ])
#
# result_irred = calculate_hueckel_secular_equation(H, info=irred, sorting_dict_values=sorting_dict_values)
#
# for i in result_irred:
#     print()
#     i.symmetry = irred
#     print(i.symmetry)
#     sp.pprint(i.eigenvector)
#     for row in range(len(i.eigenvector)):
#         print(len(i.salcs))
#         if i.eigenvector[row] != 0:
#             print({"factor": i.eigenvector[row], "salc": salcs[row]})
