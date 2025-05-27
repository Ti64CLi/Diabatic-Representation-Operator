from dataclasses import dataclass
from collections import Counter
from typing import Optional
from monome import Monome, generate_monoms
from variable import Variable
from symmetry import Symmetry

@dataclass
class ComplexInvariant:
    monome: Monome
    isreal: bool

    def is_real(self) -> bool:
        return self.is_real

    def imag(self) -> str:
        return f"ρ({str(self.monome)})"

    def real(self) -> str:
        return f"r({str(self.monome)})"

    def __str__(self) -> str:
        return f"i{self.imag()}" if self.is_real() else self.real()


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

def find_fundamental_invariants(variables: list[Variable], max_order: int, n: int) -> list[ComplexInvariant]:
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

    return [ComplexInvariant(m, m.is_real()) for m in fundamentals]

def generate_r(variables: list[Variable], n: int) -> list[ComplexInvariant]:
    return [inv for inv in find_fundamental_invariants(variables, n, n) if inv.is_real()]

def generate_irhos(variables: list[Variable], n: int) -> list[ComplexInvariant]:
    invariants = find_fundamental_invariants(variables, n, n)

    return [inv for inv in invariants if not inv.is_real()]
