from pypointgroup import PointGroup

# Punktgruppe D6h erstellen
pg = PointGroup("D6h")

# Symmetrieoperationen anzeigen
print("Symmetry operations:")
for op in pg.operations:
    print(op)

# Cayley-Tabelle anzeigen
print("\nCayley table:")
print(pg.cayley_table())
