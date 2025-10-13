

class molecule_orbital():

    def __init__(self):
        self.eigenvalue = None #Energy
        self.eigenvector = None
        self.symmetry = None
        self.occupation = None
        self.salcs = []

    def set_occupation(self, occupation):
        if occupation in [0, 1, 2]:
            self.occupation = occupation

    def get_salc_equation(self):
        eq = 0
        for x in self.salcs:
            eq += x["factor"] * x["salc"].equation
        return eq


    def getEnergy(self):
        if self.eigenvalue is not None and self.occupation is not None:
            return self.eigenvalue * self.occupation
        return "unknown"

    def print(self):
        print(f"\nMO of {self.symmetry}")
        print("- eigenvalue", self.eigenvalue)
        print("- occupation", self.occupation)
        print("- energy", self.getEnergy())
        print("- linear combination of salcs:\t", self.get_salc_equation())
        # print()


