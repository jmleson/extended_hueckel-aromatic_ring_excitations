from typing import Callable


class SymmetryOperation():
    def __init__(self, n:int, name:str, transform_p:Callable, amount:int):
        self.n = n
        self.name = name
        if amount != 1:
            raise Exception("unable to handle this later on")
        self.amount = amount

        self.transform_p = transform_p

    def transform_p(self):
        pass

    def transform_s(self):
        pass