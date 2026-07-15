from typing import Tuple
import sympy as sp

from src.main.MoleculeRepresentation import MoleculeRepresentation
from src.main.Transition import Transition, calculate_dispersion_energy
from src.main.construct_occupied_triplett_states import construct_occupied_triplett_states
from src.main.get_ground_state import get_ground_state
from src.info_document.sketch_chemical_ring import sketch_chemical_ring
from src.main.molecule_orbital import molecule_orbital
from src.main.save_latex_export import save_latex_export
# from src.main.wrappers import safe_simplify, safe_ratsimp, safe_expand


class MoleculeState:
    def __init__(self, n:int, bound_cl_to_c_positions:list[int], n_instead_of_c:list[int], heterosymbol:str = "Cl"):
        for i in bound_cl_to_c_positions:
            if i not in range(1,n+1):
                raise Exception("no valid C atom to bind Cl/other heteroatom to")
        self.n = n
        self.bound_cl_to_c_positions = bound_cl_to_c_positions
        self.heterosymbol = heterosymbol
        self.n_instead_of_c = n_instead_of_c
        self.p = MoleculeRepresentation(n=n, bound_cl_to_c_positions=bound_cl_to_c_positions, n_instead_of_c=n_instead_of_c, heterosymbol=heterosymbol)
        self.bonding_p = None

        self.s = MoleculeRepresentation(n=n, bound_cl_to_c_positions=bound_cl_to_c_positions, n_instead_of_c=n_instead_of_c, heterosymbol=heterosymbol)
        self.bonding_s = None
        self.set_up()
        # try:
        #     self.set_up()
        # except Exception as e:
        #     print(f"!!! set up failed due to {e}")

        self.p_occupation = None


    def set_up(self):
        # print("P ORBITALS")
        p_mo_orbitals = self.p.get_energy_levels()
        # for x in p_mo_orbitals:
        #     print("\t", x.symmetry, ":\t", x.eigenvalue)
        self.bonding_p = p_mo_orbitals

        # print("\nS ORBITALS")
        self.s.set_to_s_orbitals()
        s_mo_orbitals = self.s.get_energy_levels()
        # for s in s_mo_orbitals:
        #     print("\t", s.symmetry, ":\t", s.eigenvalue)
        self.antibonding_s = []
        alpha, alpha_s, beta, beta_s = sp.symbols("alpha alpha_s beta beta_s")
        for s in range(len(s_mo_orbitals)):
            s_mo_orbitals[s].eigenvalue = s_mo_orbitals[s].eigenvalue.subs(alpha, alpha_s).subs(beta, beta_s)
            if s_mo_orbitals[s].eigenvector is not None:
                s_mo_orbitals[s].eigenvector = s_mo_orbitals[s].eigenvector.subs(alpha, alpha_s).subs(beta, beta_s)
            self.antibonding_s.append(s_mo_orbitals[s])


    def get_ground_state(self):
        return get_ground_state(nel=self.n, norb=self.n) + [2]*len(self.bound_cl_to_c_positions)

    def set_occupation(self, p_occupation):
        if len(p_occupation) != self.n + len(self.bound_cl_to_c_positions):
            raise Exception("wrong number of orbitals")
        for n in p_occupation:
            if n not in [0,1,2]:
                raise Exception("wrong occupation")
        self.p_occupation = p_occupation



    def compare_ground_state_assumption(self, print_active:bool = True):
        ground_state_occupation = self.get_ground_state()
        # compare ground-state <-> ground-state assumption:
        self.set_occupation(ground_state_occupation)
        ground_state_with_warning = self.calculate_result_for_all_transitions_of_set_occupation(print_active=print_active)
        ground_state_with_warning = ground_state_with_warning[0]
        ground_state = [ground_state_with_warning[0]]
        e_disp = calculate_dispersion_energy(ground_state, ground_state)
        if print_active:
            print(f"X) A =", ground_state_occupation, "<-> B =", ground_state_occupation)
            print(set([g[1] for g in ground_state_with_warning]))
            print("\tE_dispersion =", e_disp)
        else:
            result = "Ground-State:\n\n"
            # if len(str([g[1] for g in ground_state_with_warning])) > 0:
            #     result += str([g[1] for g in ground_state_with_warning]) +"\n\n"
            result += "A =" + str(ground_state_occupation) + r"$\leftrightarrow$ B =" + str(ground_state_occupation) +"\n"
            simplified_latex_expression = save_latex_export(e_disp)
            result +=  r"$$E_{dispersion} =" +"\n"+ simplified_latex_expression + "$$\n"
            return result

    def calculate_result_for_all_transitions_results(self, print_active:bool=True):
        triplet_states = self.construct_occupied_triplet_states()
        result = "\n"+ r"\begin{enumerate}"+"\n"
        state_no = 0
        for a in triplet_states:
            for b in triplet_states:
                if a["multiplicity"] + b["multiplicity"] == 5 + 1:
                    state_no += 1
                    if print_active:
                        print(f"{state_no}) A =", a["occupation"], "<-> B =", b["occupation"])
                    else:
                        result += fr""" \item State: A = {a["occupation"]} $\leftrightarrow$ B = {b["occupation"]} """+"\n"
                    self.set_occupation(p_occupation=a["occupation"])
                    triplett_states_A_with_warning = self.calculate_result_for_all_transitions_of_set_occupation(print_active=False)
                    triplett_states_A = [t[0] for t in triplett_states_A_with_warning]
                    self.set_occupation(p_occupation=b["occupation"])
                    triplett_states_B_with_warning = self.calculate_result_for_all_transitions_of_set_occupation(print_active=False)
                    triplett_states_B = [t[0] for t in triplett_states_B_with_warning]
                    e_disp = calculate_dispersion_energy(triplett_states_A, triplett_states_B)
                    if print_active:
                        print("\tE_dispersion =", e_disp)
                    else:
                        result += r"$$E_{dispersion} =" + save_latex_export( e_disp ) + "$$\n"
                else:
                    # skip triplet quintet combinations
                    pass
                    # print("Multiplicity",a["multiplicity"] + b["multiplicity"], ", A =", a["occupation"], "<-> B=", b["occupation"] )
        result += r"\end{enumerate}"+"\n"
        return result

    def construct_occupied_triplet_states(self):
        triplet_states = construct_occupied_triplett_states(nel=self.n, norb=self.n)
        if len(self.bound_cl_to_c_positions) != 0:
            to_add = (2,) * len(self.bound_cl_to_c_positions)
            # adjust to full Cl orbital:
            for triplet in triplet_states:
                triplet["occupation"] += to_add
        return triplet_states

    def calculate_result_for_all_transitions(self, print_active:bool = True):
        triplet_states = self.construct_occupied_triplet_states()
        transitions = ""
        for triplet in triplet_states:
            if print_active:
                print(f"\033[1mState: {triplet}\033[0m")
            transitions += (r"\subsection*{State with Occupation "
                            + f"{triplet['occupation']} ({triplet['unpaired electrons']} "
                            +  f"unpaired electrons, mult {triplet['multiplicity']})"
                            + r"}")
            self.set_occupation(p_occupation=triplet["occupation"])
            allowed_transitions_with_warning = self.calculate_result_for_all_transitions_of_set_occupation(print_active=print_active)
            for i in allowed_transitions_with_warning:
                transition, warning = i
                transitions += transition.to_latex() + "\n" + warning
        return transitions

    def calculate_energy_before_transition(self):
        if self.p_occupation is None:
            raise Exception("perform set_occupation() before")
        return self.calculate_energy_for_state_and_occupation(state=self.bonding_p, occupation=self.p_occupation)

    def calculate_energy_for_state_and_occupation(self, state:list[molecule_orbital], occupation:tuple[int,...]) -> Tuple[sp.Expr, str]:
        warning = ""
        if len(state) != self.n + len(self.bound_cl_to_c_positions):
            warning = ("\t" + r"\textcolor{red}{CAUTION:} "
                    + r"state should have same length as occupation $\rightarrow$ missing results?" + "\n\t"
                    + r"DONT TRUST \glqq{}Energy of State before Excitation\grqq{}, only relative values are correct"
                    + "\n\n")
            # raise Exception("state and occupation must have same length")
        if len(state) > self.n + len(self.bound_cl_to_c_positions):# or len(occupation) != self.n + len(self.bound_cl_to_c_positions):
            raise Exception("state too long // occupation must have length of p orbitals")
        e = 0
        for i in range(len(state)):
            if len(state[i].salcs) != len(occupation):
                if len(state[i].salcs) > len(occupation):
                    raise Exception("unfitting occupation for this state!")
                if ( len(occupation) != self.n + len(self.bound_cl_to_c_positions)):
                    warning += r"less orbital occupations set than relevant for this SALC $\rightarrow$ "+ "\n"
                    warning += f"{len(occupation)} zu {self.n} + {len(self.bound_cl_to_c_positions)} \n"
                    warning += f"{len(state[i].salcs)} zu {len(occupation)}, {occupation} \n"
                    warning += (r"\textcolor{red}{DONT TRUST \glqq{}Energy of State before Excitation\grqq{}, "
                                r"only relative values are correct}")+ "\n\n"
            state[i].occupation = occupation[i]
            # state[i].print()
            e += state[i].getEnergy()
        return e, warning

    def get_p_after_transition(self, s_after_transition:tuple[int]):
        return tuple([self.p_occupation[i] - s_after_transition[i] for i in range(len(self.p_occupation))])

    def get_thinkable_transitions(self) -> Tuple[list[Transition], list[str]]:
        if self.p_occupation is None:
            raise Exception("perform set_occupation() before")
        if len(self.p_occupation) != self.n+ len(self.bound_cl_to_c_positions):
            raise Exception(f"unfitting occupation sequence for {self.n} orbitals")
        warnings = []
        e_state_before, warning = self.calculate_energy_before_transition()
        warnings.append(warning)
        # 1st orbital (lowest energy) unchanged
        s_occupations = []
        for i in range(1, len(self.p_occupation)):
            if self.p_occupation[i] > 0 and i in range(self.n):
                s_after_transition = tuple([ 1 if p == i else 0 for p in range(len(self.p_occupation)) ])
                transition = Transition(n=self.n+len(self.bound_cl_to_c_positions),
                                        s_after_transition=s_after_transition,
                                        p_before_transition=self.p_occupation)

                E_s_state, warning_state1 = self.calculate_energy_for_state_and_occupation(state=self.antibonding_s,
                                                                   occupation=s_after_transition)
                energy_single_p, warning_state2 = self.calculate_energy_for_state_and_occupation(state=self.bonding_p,
                                                                   occupation=s_after_transition)
                warnings.append(warning_state1+" - "+warning_state2)
                E_state_after = E_s_state + (e_state_before - energy_single_p)
                changed_orbital_index = transition.get_changed_orbital_index()
                try:
                # if True:
                    energy_zero_p, warning_zero_p = self.calculate_mean_energy_for_state(state=self.bonding_p)
                    energy_zero_s, warning_zero_s = self.calculate_mean_energy_for_state(state=self.antibonding_s)
                    warnings[-1] += warning_zero_p + " AND " + warning_zero_s + "\n"
                    transition.set_up(orbital_to_excite_of = self.bonding_p[changed_orbital_index],
                                  orbital_to_excite_to = self.antibonding_s[changed_orbital_index],
                                  energy_zero_p = energy_zero_p,
                                  energy_zero_s = energy_zero_s,
                                  energy_of_state_before_excitation = e_state_before,
                                  energy_of_state_after_excitation = E_state_after)

                    assert E_s_state - energy_single_p == E_state_after - e_state_before
                    assert transition.get_transitioning_energy() == E_state_after - e_state_before

                    s_occupations.append(transition)
                except IndexError as e:
                    print("\t something flaut with transition",
                          "(probably due to unsolved H-matrices and thereby missing p orbital energies)",
                          changed_orbital_index, ":", e)
        return s_occupations, warnings

    def calculate_mean_energy_for_state(self, state) -> Tuple[sp.Expr, str]:
        # calculate "zero" energy of orbital set, as defined in state (p_bonding / s_antibonding)
        energy_of_single_total_occupation, warning = self.calculate_energy_for_state_and_occupation(
                                                        state=state,
                                                        occupation=tuple([1 for i in range(len(state))])
                                                    )
        return energy_of_single_total_occupation / len(state), warning

    def symmetry_allowed_transition(self, transition:Transition):
        sym_p = transition.orbital_to_excite_of.symmetry
        sym_s = transition.orbital_to_excite_to.symmetry
        dipole_operator = self.p.pointgroup.dipole_operator_symmetry
        # calculate symmetry of total integral: sym_p  x  dipole_operator  x  sym_s
        part_I = self.p.pointgroup.multiply(dipole_operator, sym_s)# new pseudo p orbital
        part_II = [self.p.pointgroup.multiply(sym_p, i) for i in part_I]
        # print(dipole_operator, "x", sym_s, "x", sym_p,  "=", part_II)
        return self.p.pointgroup.total_symmetric_representation in [x for row in part_II for x in row]

    def calculate_result_for_all_transitions_of_set_occupation(self, print_active:bool=True):
        result = self.get_thinkable_transitions()
        transitions, warnings = result

        allowed_transitions = []
        for transition, warning in zip(transitions, warnings):
            if self.symmetry_allowed_transition(transition):
                if print_active:
                    transition.print()
                    print(warning)
                allowed_transitions.append((transition, warning))
            else:
                if print_active:
                    print("\tFORBIDDEN\n")
            #     transition.print()
        return allowed_transitions



    def sketch_chemical_ring(self, filename = None):
        atom_symbols = [
            "N" if i in self.n_instead_of_c else "C"
            for i in range(self.n)
        ]
        if filename is None:
            filename = f"molecule_{self.n}.png"
        return sketch_chemical_ring(speichername=filename,
                             atom_symbols=atom_symbols, heterosymbol = self.heterosymbol,
                             bound_Cl_to_C=self.bound_cl_to_c_positions)

    def latex_datei_erstellen(self, molekuel_name):
        """Erstellt eine einfache LaTeX-Datei mit Molekülname, Bild und Punktgruppe."""
        tex_datei = molekuel_name.replace(" ","").lower()+".tex"
        header = (fr"""
        \documentclass[12pt,a4paper]{{article}}
        \usepackage{{graphicx}}
        \usepackage{{geometry}}
        \geometry{{margin=2cm}}
        \usepackage{{helvet}}
        \usepackage{{xcolor}}
        \usepackage{{amsmath}} % for bmatrix
        \setlength{{\parindent}}{{0pt}}  % No indentation globally
        \renewcommand{{\familydefault}}{{\sfdefault}}
        \usepackage[german]{{babel}}% for glqq / grqq

        \begin{{document}}
        """.strip())
        content = "" #INFO added later on
        footer = "\n\n"+fr"""
        \end{{document}}
        """.strip()

        try:
            bild_datei = self.sketch_chemical_ring(filename = molekuel_name.replace(" ","").lower())
            punktgruppe = self.p.pointgroup.__class__.__name__
            content += fr"""
            \begin{{center}}
                \Huge \textbf{{{molekuel_name}}} \\[1cm]
                \includegraphics[width=0.5\textwidth]{{{bild_datei}}} \\[0.5cm]
                \Large Point Group: \textbf{{{punktgruppe}}}
            \end{{center}}
            """.strip()
        except Exception as e:
            raise Exception(e)
        content += r"\newpage"

        try:
            content += r"\section*{p orbitals}"
            content += self.p.get_latex_symmetry_behavior()+"\n\n"
            content += self.p.get_latex_salcs(print_active=False)+"\n\n"

            content += r"\newpage \section*{s orbitals}"
            content += self.s.get_latex_symmetry_behavior() + "\n\n"
            content += self.s.get_latex_salcs(print_active=False) + "\n\n"
        except Exception as e:
            msg = str(e).replace("_", r"\_")
            content += r"\textcolor{red}{"+f"Error in solving salc-hueckel-matrix {msg}\n"+"}"

        content += r"\newpage \section*{Transitions from p orbitals into s orbitals}"
        try:
            content += self.calculate_result_for_all_transitions(print_active=False)
        except Exception as e:
            content += r"\textcolor{red}{"+f"unable to calculate transitions due to: {e}"+"}"

        content += r"\newpage \section*{Dispersion energies following from the transitions}"+ "\n"
        try:
            content += self.calculate_result_for_all_transitions_results(print_active=False)
        except Exception as e:
            msg = str(e).replace("_", r"\_")
            content += r"\textcolor{red}{"+f"Error in Calculations Dispersion {msg}"+"}"
        try:
            x = self.compare_ground_state_assumption(print_active=False)
            content += x
        except Exception as e:
            msg = str(e).replace("_", r"\_")
            content += r"\textcolor{red}{"+f"Error for Ground State {msg}"+"}"

        content = content.replace(r"\alpha", r"\tilde{\alpha}").replace(r"\beta", r"\tilde{\beta}").replace(r"\delta",r"\tilde{\delta}")
        content = content.replace(r"\Gamma", r"\mathsf{\Gamma}")
        content = content.replace(r"p_", r"\text{p}_").replace(r"s_", r"\text{s}_")
        with open("RESULTS/"+tex_datei, "w") as f:
            f.write(header + content + footer)
        print(f"LaTeX-Datei gespeichert als: {tex_datei}")
        return
