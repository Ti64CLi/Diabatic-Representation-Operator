from collections import Counter
from itertools import combinations_with_replacement
from monome import Monome
from variable import Variable
from invariant import ComplexInvariant

def filter_appearing_variables(variables: list[Variable]) -> list[Variable]:
    """
    Filter variables to keep only appearing ones (excludes A1 variables)
    """
    avariables = []

    for v in variables:
        if not v.symmetry.is_A1():
            avariables.append(v)

    return avariables

def try_to_factorize(monome: Monome, factors: list[Monome]) -> bool:
    target = Counter(monome.variables)
    degree = target.total()

    for v in target:
        if v.symmetry.is_A2() and target[v] >= 2:
            return True

    for factor in factors:
        combo = Counter(factor.variables)
        combo_degree = combo.total()

        if degree < combo_degree:
            continue

        if combo <= target:
            return True

    return False

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

        if monoms[m] != 0:
            m.complex_conjugate = False
        else:
            if monoms[m.conjugate()] != 0:
                m.complex_conjugate = True
            else:
                m.complex_conjugate = False


        if remove_cc and m.complex_conjugate == True:
            m = m.conjugate()

        if monoms[m] == 0:
            monoms[m] = 1

    return list(monoms.elements())

def find_fundamental_invariants(variables: list[Variable], max_order: int, n: int, remove_cc: bool = True) -> list[ComplexInvariant]:
    """
    Find the fundamental invariants, those which cannot be factorize anymore
    """
    # avariables = filter_appearing_variables(variables)
    fundamentals = []

    for order in range(1, max_order + 1):
        monoms = generate_monoms(variables, order, n, remove_factorizable=False, remove_cc=remove_cc)

        for m in monoms:
            if m.weight_mod(n) != 0:
                continue

            if not try_to_factorize(m, fundamentals):
                fundamentals.append(m)

    return [ComplexInvariant(finv, finv.is_real()) for finv in fundamentals]

def generate_appearing_monoms(variables: list[Variable], n: int, remove_cc: bool = True) -> list[Monome]:
    """
    Generate all appearing monomials
    """
    avariables = filter_appearing_variables(variables)
    finvs = list(map(lambda x:x.monome, find_fundamental_invariants(avariables, n, n, remove_cc=remove_cc)))
    monoms = []
    amonoms = []

    for order in range(n):
        monoms.extend(generate_monoms(avariables, order, n, remove_factorizable=True, remove_cc=remove_cc))

    for m in monoms:
        if not try_to_factorize(m, finvs):
            amonoms.append(m)

    return amonoms
