import matplotlib.pyplot as plt
import numpy as np

def sketch_chemical_ring(speichername=None,
                         atom_symbols=None, bound_Cl_to_C=None):
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
    radius = 1
    abstand_Cl = 0.5
    farbe = 'blue'#'darkorange'
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
    ax.plot(x_plot, y_plot, color=farbe, linewidth=2)

    # Atomsymbole einzeichnen
    for i, (xi, yi, atom_symbol) in enumerate(zip(x, y, atom_symbols)):
        ax.text(xi, yi, atom_symbol, fontsize=14, fontweight='bold',
                ha='center', va='center')

        # Falls Cl an dieses Atom gebunden ist
        if bound_Cl_to_C and i in bound_Cl_to_C:
            # Richtung nach außen (vom Mittelpunkt weg)
            richtung = np.array([xi, yi]) / np.linalg.norm([xi, yi])
            Cl_pos = np.array([xi, yi]) + richtung * abstand_Cl
            ax.plot([xi, Cl_pos[0]], [yi, Cl_pos[1]], color=farbe, linewidth=1.5)
            ax.text(Cl_pos[0], Cl_pos[1], 'Cl', fontsize=12, fontweight='bold',
                    ha='center', va='center', color='green')

    ax.set_aspect('equal', 'box')
    ax.axis('off')

    if speichername:
        plt.savefig(speichername, bbox_inches='tight', transparent=True)
        print(f"Bild gespeichert als: {speichername}")
    else:
        plt.show()
    return speichername


if __name__ == '__main__':
    sketch_chemical_ring(
        atom_symbols=['C'] * 3,
        bound_Cl_to_C=[1, 4],
        speichername='hexagon.png'
    )
