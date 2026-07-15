import matplotlib.pyplot as plt
import numpy as np

def sketch_chemical_ring(speichername=None,
                         atom_symbols=None, bound_Cl_to_C=None, heterosymbol="Cl"):
    """
    Zeichnet ein n-Eck mit angegebenen Atomsymbolen und optional angehängten Cl-Atomen.

    Parameters
    ----------
    radius : float
        Radius des Rings.
    farbe : str
        Farbe der Linien.
    speichername : str oder None
        Wenn angegeben, wird das Bild gespeichert.
    atom_symbols : list[str]
        Liste der Atomsymbole, z. B. ['C', 'C', 'C', 'C', 'C', 'C'].
    bound_Cl_to_C : list[int] oder None
        Indizes (0-basiert), an denen Cl-Atome an die angegebenen C gebunden werden sollen.
    abstand_Cl : float
        Abstand der Cl-Atome nach außen.
    """
    distances = get_p_distances(speichername)
    box = dict(
                    boxstyle="round,pad=0.1",
                    facecolor="white",
                    edgecolor="none"
                )

    radius = 1
    abstand_Cl = 0.7
    farbe = 'grey'#'darkorange'
    if atom_symbols is None:
        atom_symbols = ['C'] * 6  # Default: Benzol
    n = len(atom_symbols)
    winkel = np.linspace(0, 2 * np.pi, n, endpoint=False)

    x = radius * np.cos(winkel)
    y = radius * np.sin(winkel)

    # Figur schließen
    x_plot = np.append(x, x[0])
    y_plot = np.append(y, y[0])

    fig, ax = plt.subplots()
    ax.plot(x_plot, y_plot, color=farbe, linewidth=5)

    # Atomsymbole einzeichnen
    for i, (xi, yi, atom_symbol) in enumerate(zip(x, y, atom_symbols)):
        ax.text(xi, yi, f"{atom_symbol}", fontsize=30, fontweight='bold',
                ha='center', va='center')
        ax.text(xi+distances["C_x"], yi-distances["C_y"], fr"($\mathbf{{p_{i+1}}}$)", fontsize=20, fontweight='bold',
                ha='center', va='center')


        # Falls Cl an dieses Atom gebunden ist
        if bound_Cl_to_C and i+1 in bound_Cl_to_C:
            # Richtung nach außen (vom Mittelpunkt weg)
            richtung = np.array([xi, yi]) / np.linalg.norm([xi, yi])
            Cl_pos = np.array([xi, yi]) + richtung * abstand_Cl
            ax.plot([xi, Cl_pos[0]], [yi, Cl_pos[1]], color=farbe, linewidth=5)
            ax.text(Cl_pos[0], Cl_pos[1], s=heterosymbol, fontsize=30, fontweight='bold',
                    ha='center', va='center', color='green')
            ax.text(Cl_pos[0]+distances["Hetero_x"], Cl_pos[1]-distances["Hetero_y"], fr"($\mathbf{{p_{{{i+1+len(atom_symbols)}}}}}$)", fontsize=20, fontweight='bold',
                ha='center', va='center', color='green',
                    )

    ax.set_aspect('equal', 'box')
    ax.axis('off')

    if speichername:
        plt.savefig("RESULTS/"+speichername, bbox_inches='tight', transparent=True, pad_inches=0.5)
        print(f"Picture saved: RESULTS/{speichername}")
    else:
        plt.show()
    return speichername



def get_p_distances(speichername:str):
    speichername = speichername.replace(".png","").replace(".pdf", "").lower()

    distances = {"C_x": 0, "C_y": 0, "Hetero_x": 0, "Hetero_y": 0 }
    if speichername == "benzene":
        distances["C_x"] = 0.215
        distances["C_y"] = 0.1
    elif "chlorobenzene" in speichername and not "hexa" in speichername:
        distances["C_x"] = 0.26
        distances["C_y"] = 0.1
        distances["Hetero_x"] = 0.3
        distances["Hetero_y"] = 0.1
    elif "hexa" in speichername:
        distances["C_x"] = 0.36
        distances["C_y"] = 0.1
        distances["Hetero_x"] = 0.37#0.4
        distances["Hetero_y"] = 0.2

    return distances


if __name__ == '__main__':

    sketch_chemical_ring(
        atom_symbols=['C'] * 6,
        bound_Cl_to_C=[],
        speichername='benzene.png',
    )

    sketch_chemical_ring(
        atom_symbols=['C'] * 6,
        bound_Cl_to_C=[1],
        speichername='chlorobenzene.png',
    )

    sketch_chemical_ring(
        atom_symbols=['C'] * 6,
        bound_Cl_to_C=[0, 2, 1, 3, 4, 5, 6],
        speichername='hexachlorobenzene.png',
    )

    sketch_chemical_ring(
        atom_symbols=['C'] * 6,
        bound_Cl_to_C=[0, 2, 1, 3, 4, 5, 6],
        speichername='hexafluorobenzene.png',
        heterosymbol="F"
    )


    sketch_chemical_ring(
        atom_symbols=['C'] * 6,
        bound_Cl_to_C=[1],
        speichername='chlorobenzene-hueckelorientation.pdf',
    )
