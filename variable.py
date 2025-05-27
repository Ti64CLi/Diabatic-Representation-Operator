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
        return Variable(self.name, self.symmetry, complex_conjugate=not self.complex_conjugate)

    def weight(self) -> int:
        sign = -1 if self.complex_conjugate else 1

        return sign * self.symmetry.weight()
