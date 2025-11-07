# from MoleculeState import MoleculeState
#
#
# def latex_datei_erstellen(molekuel_name, bild_datei, punktgruppe, tex_datei):
#     """Erstellt eine einfache LaTeX-Datei mit Molekülname, Bild und Punktgruppe."""
#     header = fr"""
#     \documentclass[12pt,a4paper]{{article}}
#     \usepackage{{graphicx}}
#     \usepackage{{geometry}}
#     \geometry{{margin=2cm}}
#     \usepackage{{helvet}}
#     \renewcommand{{\familydefault}}{{\sfdefault}}
#
#     \begin{{document}}
#
#     \begin{{center}}
#         \Huge \textbf{{{molekuel_name}}} \\[1cm]
#         \includegraphics[width=0.5\textwidth]{{{bild_datei}}} \\[0.5cm]
#         \Large Punktgruppe: \textbf{{{punktgruppe}}}
#     \end{{center}}
#
#     \end{{document}}
#     """.strip()
#
#     footer = fr"""
#     \end{{document}}
#     """.strip()
#
#     with open(tex_datei, "w") as f:
#         f.write(content)
#
#     print(f"LaTeX-Datei gespeichert als: {tex_datei}")
#
#
#
#
# if __name__ == '__main__':
#     m = MoleculeState(n=6, bound_cl_to_c_positions=[], n_instead_of_c=[])
#     m.latex_datei_erstellen("Benzene")
#     # triplet_states = construct_occupied_triplett_states(nel=6, norb=6)
#     # calculate_result_for_all_transitions_print(triplet_states=triplet_states, m=m)
#     # calculate_result_for_all_transitions_results(m=m, triplet_states=triplet_states)
#     #
#     # ground_state_occupation = get_ground_state(nel=6, norb=6)
#     # compare_ground_state_assumption(m=m, ground_state_occupation=ground_state_occupation)
#
#     # m.sketch_chemical_ring()
#     #
#     # latex_datei_erstellen(
#     #     molekuel_name="Benzen",
#     #     bild_datei=f"molecule_{m.n}.png",
#     #     punktgruppe=m.p.pointgroup.__class__.__name__,
#     #     tex_datei="benzol.tex"
#     # )