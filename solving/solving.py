# import sympy as sp
#
# from solving.solve_for_eigenvalues_first import solve_for_eigenvalues_first
#
#
# def solving(H: sp.Matrix, info:str):
#     eigen_data = None  # calculate_with_timeout(H.eigenvects, (), 600, msg="Calculation of H.eigenvects() timed out.")
#     if eigen_data is None:
#         eigen_data = solve_for_eigenvalues_first(H=H)
#     # if eigen_data is None:
#     #     if H.rows == H.cols and H.rows == 2:
#     #             eigen_data = solve_2x2(H=H)
#     # if eigen_data is None:
#     #     eigen_data = alternative_solver(H=H, info=info)
#     print("Eigen data (", info, "):", "\n\t", type(eigen_data), flush=True)
#     return eigen_data