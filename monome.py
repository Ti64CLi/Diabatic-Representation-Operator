from dataclasses import dataclass
from collections import Counter
from itertools import combinations_with_replacement
from variable import Variable
from symmetry import Symmetry
from utils import *

@dataclass
class Monome:
    variables: list[Variable]

    def weight(self) -> int:
        return sum(v.weight() for v in self.variables)

    def weight_mod(self, n: int) -> int:
        return self.weight() % n

    def is_factorisable(self, n: int) -> bool:
        return len(self.variables) > 0 and self.weight_mod(n) == 0

    def is_real(self):
        if len(self.variables) == 0:
            return True

        cv = Counter(self.variables)

        for v in cv:
            if cv[v] != cv[v.conjugate()]:
                return False

        return True

    def __str__(self) -> str:
        cv = Counter(self.variables)

        if len(+cv) == 0:
            return "1"

        return ''.join((f"{v}{num2sup(cv[v])}" if cv[v] > 1 else f"{v}") for v in cv)

def generate_monoms(variables: list[Variable], order: int, n: int, remove_factorizable: bool = True) -> list[Monome]:
    """
    Generates monoms of given order and filters factorizable monoms if needed
    """
    monoms = []

    for combo in combinations_with_replacement(variables, order):
        m = Monome(combo)

        if remove_factorizable and m.is_factorisable(n):
            # add treatement for invariant factors ?
            continue

        monoms.append(m)

    return monoms

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

def try_to_factorize(monome: Monome, factors: list[Monome]) -> bool:
    target = Counter(monome.variables)
    degree = target.total()

    for factor in factors:
        combo = Counter(factor.variables)
        combo_degree = combo.total()

        if degree < combo_degree:
            continue

        if combo <= target:
            return True

    return False

def find_fundamental_invariants(variables: list[Variable], max_order: int, n: int) -> list[Monome]:
    """
    Find the fundamental invariants, those which cannot be factorize anymore
    """
    fundamentals = []

    for order in range(1, max_order + 1):
        monoms = generate_monoms(variables, order, n, remove_factorizable=False)

        for m in monoms:
            if m.weight_mod(n) != 0:
                continue

            if not try_to_factorize(m, fundamentals):
                fundamentals.append(m)

    return [Monome(variables=())] + fundamentals

def generate_appearing_monoms(variables: list[Variable], n: int) -> list[Monome]:
    """
    Generate all appearing monomials
    """
    monoms = []

    for order in range(n):
        monoms.extend(generate_monoms(variables, order, n, True))

    return monoms

def generate_rhos(variables: list[Variable], n: int) -> list[Monome]:
    invariants = find_fundamental_invariants(variables, n, n)

    return [Monome(variables=())] + [inv for inv in invariants if not inv.is_real()]
