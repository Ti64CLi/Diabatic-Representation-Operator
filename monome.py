from dataclasses import dataclass
from collections import Counter
from itertools import combinations_with_replacement
from variable import Variable
from symmetry import Symmetry
from utils import *
from invariant import InvariantType

@dataclass
class Monome:
    variables: list[Variable]
    complex_conjugate: bool = False
    invariant_type: InvariantType = None

    def __str__(self) -> str:
        cv = +Counter(self.variables)

        if len(cv) == 0: # shouldn't have any negative
            return "1"
        else:
            s = '*'.join((f"{v}{num2sup(cv[v])}" if cv[v] > 1 else f"{v}") for v in cv)

            if self.invariant_type is not None:
                if self.invariant_type.is_invariant():
                    if self.invariant_type.is_real_invariant():
                        return f"r({s})"
                    elif self.invariant_type.is_pseudo_invariant():
                        return f"ρ²({s})"
                    else: # full invariant
                        return s
                elif self.invariant_type.imag == True:
                    return f"iρ({s})"
            else:
                return s

        return ""

    def __eq__(self, other) -> bool:
        assert isinstance(other, Monome)

        return Counter(self.variables) == Counter(other.variables)
        # should we test for equality of complex_conjugate ?
        # it's supposed to be equal if variables are equal

    def __hash__(self) -> int:
        # doesn't integrate self.complex_conjugate because of its variant character
        return hash(str(self.variables) + str(self.invariant_type))

    def weight(self) -> int:
        return sum(v.weight() for v in self.variables)

    def weight_mod(self, n: int) -> int:
        return self.weight() % n

    def is_Cn_invariant(self, n: int) -> bool:
        return self.weight_mod(n) == 0

    def is_sigman_invariant(self, n: int) -> bool:
        return Counter(self.variables) == Counter(self.conjugate().variables)

    def is_factorisable(self, n: int) -> bool:
        ab1_var = sum(1 for v in self.variables if v.symmetry.is_A1() or v.symmetry.is_B1())
        a2_var = sum(1 for v in self.variables if v.symmetry.is_A2())
        b2_var = sum(1 for v in self.variables if v.symmetry.is_B2())

        return len(self.variables) == 0 or (self.weight_mod(n) == 0 and a2_var != 1 and b2_var != 1) or self.weight() > n or ab1_var > 0

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

    def is_pure_imag(self):
        if len(self.variables) == 0:
            return False

        cv = Counter(self.variables)
        ab2_sym = 0

        for v in cv:
            if v.symmetry.is_A2() or v.symmetry.is_B2():
                ab2_sym += 1

        return ab2_sym % 2 == 1

    def conjugate(self):
        ccvariables = []

        for v in self.variables:
            ccvariables.append(v.conjugate())

        return Monome(ccvariables, not self.complex_conjugate, self.invariant_type)
