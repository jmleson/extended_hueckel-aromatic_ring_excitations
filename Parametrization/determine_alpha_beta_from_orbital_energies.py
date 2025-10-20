
import numpy as np

def determine_alpha_beta_from_orbital_energies(orbital_energies_increasing_E:list[float]):
    alpha = np.mean(orbital_energies_increasing_E)
    print()
    print(f"alpha = {alpha:.6f}")

    beta_6 = (alpha - orbital_energies_increasing_E[6 - 1] )/2
    beta_5 = alpha - orbital_energies_increasing_E[5 - 1]
    beta_4 = alpha - orbital_energies_increasing_E[4 - 1]
    beta_3 = -alpha + orbital_energies_increasing_E[3 - 1]
    beta_2 = -alpha + orbital_energies_increasing_E[2 - 1]
    beta_1 = (-alpha + orbital_energies_increasing_E[1 - 1]) / 2
    beta_proposals = [beta_1, beta_2, beta_3, beta_4, beta_5, beta_6]
    # print(beta_proposals)
    beta = np.mean(beta_proposals)
    print(f"beta = {beta:.6f}")
    print()
    return alpha, beta





if __name__ == "__main__":
    # Benzene, RHF-orbitals, D2h:
    phi_1 = -0.50520# 1.5
    phi_2 = -0.33696# 1.7
    phi_3 = -0.33696# 1.6
    phi_4 = 0.13153# 2.5
    phi_5 = 0.13153# 1.8
    phi_6 = 0.33593# 2.7
    alpha, beta = determine_alpha_beta_from_orbital_energies([phi_1, phi_2, phi_3, phi_4, phi_5, phi_6])



    ''' FRAGLICHE AUSWAHL!!! 
    10.1: {'7', '5', '1', '11'}, 0.41664 	 {'s': -0.5165800000000003, 'px': 0.45433, 'py': -0.10866000000000003} 

	 8.1: {'7', '5', '1', '11'}, 0.20244 	 {'s': 0.23784000000000027, 'px': 1.20323, 'py': 4.20892} 

	 4.4: {'7', '5', '1'}, 0.20244 	 {'s': 1.9348899999999998, 'py': -1.20324, 'px': -3.64189} 

	 5.2: {'7', '5', '1'}, 0.17014 	 {'s': 0.5855000000000001, 'px': -1.3728600000000002, 'py': -1.41959} 

	 6.3: {'7', '5', '1', '11'}, 0.17013 	 {'s': 0.1400699999999997, 'py': 1.7138399999999998, 'px': 1.41958} 

	 7.1: {'7', '5', '1', '11'}, 0.14657 	 {'s': 0.6292399999999998, 'px': -0.56874, 'py': 0.13602000000000003} 
	 '''
    xi_1 = 0.41664
    xi_2 = 0.20244
    xi_3 = 0.20244
    xi_4 = 0.17014
    xi_5 = 0.17013
    xi_6 = 0.14657
    alpha, beta = determine_alpha_beta_from_orbital_energies([xi_1, xi_2, xi_3, xi_4, xi_5, xi_6])