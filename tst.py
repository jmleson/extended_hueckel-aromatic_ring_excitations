from molsym import PointGroup

# Punktgruppe D6h erstellen
pg = PointGroup("D6h")

# Alle irreduziblen Darstellungen abfragen
for irrep in pg.irreps:
    print(f"Name: {irrep.name}")
    print(f"Dimension: {irrep.dim}")
    print(f"Matrixdarstellungen: {irrep.mats}")
    print()
