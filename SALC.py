from IrreducibleRepresentation import IrreducibleRepresentation


class SALC:
    def __init__(self, n:int, irred:IrreducibleRepresentation):
        self.n = n
        self.irred = irred
        self.prefactors_of_AOs = []

