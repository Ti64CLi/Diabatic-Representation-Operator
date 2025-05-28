from dataclasses import dataclass
from collections import Counter
from itertools import combinations_with_replacement
from variable import Variable
from symmetry import Symmetry
from utils import *

@dataclass
class Monome:
    variables: list[Variable]
    complex_conjugate: bool = False

    def __str__(self) -> str:
        cv = Counter(self.variables)

        if len(+cv) == 0: # shouldn't have any negative
            return "1"

        return ''.join((f"{v}{num2sup(cv[v])}" if cv[v] > 1 else f"{v}") for v in cv)

    def __eq__(self, other) -> bool:
        assert isinstance(other, Monome)

        return Counter(self.variables) == Counter(other.variables)
        # should we test for equality of complex_conjugate ?
        # it's supposed to be equal if variables are equal

    def __hash__(self):
        # doesn't integrate self.complex_conjugate because of its invariant character
        return hash(str(self.variables))

    def weight(self) -> int:
        return sum(v.weight() for v in self.variables)

    def weight_mod(self, n: int) -> int:
        return self.weight() % n

    def is_factorisable(self, n: int) -> bool:
        return len(self.variables) > 0 and (self.weight_mod(n) == 0 or self.weight() > n)

    def is_real(self):
        if len(self.variables) == 0:
            return True

        cv = Counter(self.variables)
        a2_sym = 0
        b2_sym = 0

        for v in cv:
            if v.symmetry.is_E() and cv[v] != cv[v.conjugate()]:
                return False
            if v.symmetry.is_A2():
                a2_sym += 1
            if v.symmetry.is_B2():
                b2_sym += 1

        return True and (a2_sym % 2 == 0) and (b2_sym % 2 == 0)

    def conjugate(self):
        ccvariables = []

        for v in self.variables:
            ccvariables.append(v.conjugate())

        return Monome(ccvariables, not self.complex_conjugate)
