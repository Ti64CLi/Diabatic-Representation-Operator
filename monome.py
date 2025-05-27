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
        return hash(str(self.variables) + str(self.complex_conjugate))

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

        for v in cv:
            if cv[v] != cv[v.conjugate()]:
                return False

        return True

    def conjugate(self):
        ccvariables = []

        for v in self.variables:
            ccvariables.append(v.conjugate())

        return Monome(ccvariables, not self.complex_conjugate)


def generate_monoms(variables: list[Variable], order: int, n: int, remove_factorizable: bool = True, remove_cc: bool = True) -> list[Monome]:
    """
    Generates monoms of given order and filters factorizable monoms if needed
    """
    monoms = Counter()

    for combo in combinations_with_replacement(variables, order):
        m = Monome(list(combo))

        if remove_factorizable and m.is_factorisable(n):
            # add treatement for invariant factors ?
            continue

        if m.weight() < 0:
            m.complex_conjugate = True
        else:
            m.complex_conjugate = False

        if remove_cc and m.complex_conjugate == True:
            m = m.conjugate()

        if monoms[m] == 0:
            monoms[m] = 1

    return list(monoms.elements())

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

def generate_appearing_monoms(variables: list[Variable], n: int, remove_cc: bool = True) -> list[Monome]:
    """
    Generate all appearing monomials
    """
    avariables = []

    for v in variables:
        if v.symmetry.is_B() or v.symmetry.is_E():
            avariables.append(v)

    monoms = []

    for order in range(n):
        monoms.extend(generate_monoms(avariables, order, n, remove_factorizable=True, remove_cc=remove_cc))

    return monoms
