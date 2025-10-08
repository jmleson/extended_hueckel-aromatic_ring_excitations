

class molecule_orbital():

    def __init__(self):
        self.eigenvalue = 0 #Energy
        self.eigenvector = None
        self.symmetry = None
        self.occupation = 0

    def set_occupation(self, occupation):
        if occupation in [0, 1, 2]:
            self.occupation = occupation


