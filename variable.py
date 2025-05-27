from dataclasses import dataclass
from symmetry import Symmetry
from utils import *

@dataclass(frozen=True)
class Variable:
    name: str
    symmetry: Symmetry
    complex_conjugate: bool = False

    def __str__(self) -> str:
        s = f"{self.name}{sign2sub(-1 if self.complex_conjugate else 1)}"

        return s + (num2sub(self.symmetry.gamma) if self.symmetry.is_E() and self.symmetry.gamma > 1 else "")

    def __eq__(self, other) -> bool:
        assert isinstance(other, Variable)

        return self.name == other.name and self.symmetry == other.symmetry # is complex_conjugate equal or not ?

    def conjugate(self):
        if not self.symmetry.is_E():
            return self

        return Variable(self.name, self.symmetry, complex_conjugate=not self.complex_conjugate)

    def weight(self) -> int:
        sign = -1 if self.complex_conjugate else 1

        return sign * self.symmetry.weight()

    def is_real(self) -> bool:
        return self.symmetry.is_A1() or self.symmetry.is_B1()

def generate_variables_list(nvarsym: list[int], n: int) -> list[Variable]:
    """
    Generates a list of variables given the number of variables of each symmetry
    """
    if len(nvarsym) == 0:
        return []

    variables = []

    # generate A1 variables
    for i in range(nvarsym[0]):
        variables.append(Variable(f"R{i + 1}", Symmetry("A1")))

    # generate A2 variables
    if len(nvarsym) > 1:
        for i in range(nvarsym[1]):
            variables.append(Variable(f"ρ{i + 1}", Symmetry("A2")))

    # generate B1 variables
    if len(nvarsym) > 2:
        if nvarsym[2] > 0:
            assert n % 2 == 0

        for i in range(nvarsym[2]):
            variables.append(Variable(f"Rb{i + 1}", Symmetry("B1", gamma=n//2)))

    # generate B2 variables
    if len(nvarsym) > 3:
        if nvarsym[3] > 0:
            assert n % 2 == 0

        for i in range(nvarsym[3]):
            variables.append(Variable(f"ρb{i + 1}", Symmetry("B2", gamma=n//2)))

    # generate E variables
    if len(nvarsym) > 4:
        for gamma3 in range(4, len(nvarsym)):
            for i in range(nvarsym[gamma3]):
                name = f"Q{i + 1}"
                sym = Symmetry("E", gamma=gamma3 - 3)
                variables.append(Variable(name, sym))
                variables.append(Variable(name, sym, complex_conjugate=True)) # add Q+ and Q-

    return variables
