
import sympy as sp

from calculate_hueckel_secular_equation import calculate_hueckel_secular_equation

# Define the symbols
E, alpha, alpha_N, beta, beta_N = sp.symbols('E alpha alpha_N beta beta_N')

# Define the matrix
matrix =  sp.Matrix([
            [alpha_N,               beta_N * sp.sqrt(2),    0,                  0],
            [beta_N * sp.sqrt(2),   alpha,                  beta,               0],
            [0,                     beta ,                  alpha,              sp.sqrt(2)*beta ],
            [0,                     0,                      sp.sqrt(2)*beta,    alpha ]
        ])

calculate_hueckel_secular_equation(H=matrix,info="bla",sorting_dict_values={})



