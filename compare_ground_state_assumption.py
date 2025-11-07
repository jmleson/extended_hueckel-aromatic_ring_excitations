from MoleculeState import MoleculeState
from Transition import calculate_dispersion_energy


# def compare_ground_state_assumption(m:MoleculeState, ground_state_occupation:list[int]):
#     # compare ground-state <-> ground-state assumption:
#     m.set_occupation(ground_state_occupation)
#     ground_state = m.calculate_result_for_all_transitions_of_set_occupation(print_active=False)
#     e_disp = calculate_dispersion_energy(ground_state, ground_state)
#     print(f"X) A =", ground_state_occupation, "<-> B =", ground_state_occupation)
#
#     # import sympy as sp
#     # alpha, beta, alpha_s, beta_s, alpha_Cl, beta_Cl, alpha_s_Cl, beta_s_Cl = sp.symbols(
#     #     f"alpha beta alpha_s beta_s alpha_Cl beta_Cl alpha_s_Cl beta_s_Cl")
#     #
#     # e_disp = e_disp.subs(-2 * alpha + 2 * alpha_s, sp.Symbol(f"B"))
#     # e_disp = e_disp.subs(- alpha + alpha_s - beta + beta_s, sp.Symbol(f"A"))
#     print("\tE_dispersion =", e_disp.simplify())
