from dataclasses import dataclass
from collections import Counter
from itertools import combinations_with_replacement
from monome import Monome
from variable import Variable
from invariant import InvariantType
from utils import sign, num2sup

type ExpansionTerm = dict[Monome, complex]

@dataclass
class MonomialExpansion:
    # [order -> [Monome -> count]]
    # ex : [1 -> [Q+ -> +i, Q+2 -> -1]]
    expansion: dict[int, ExpansionTerm]

    def __str__(self) -> str:
        if len(self.expansion) == 0:
            return "0"

        s = ""

        for order, exp in self.expansion.items():
            if order == 0 and len(exp) != 0:
                s += "+1"

                continue

            for monome, coeff in exp.items():
                if coeff != 0:
                    if coeff.real != 0:
                        s += sign(coeff.real)
                        s += f"Re(({monome}){num2sup(order)})"

                    if coeff.imag != 0:
                        s += sign(coeff.imag)
                        s += f"Im(({monome}){num2sup(order)})"

        return s

    def __add__(self, other):
        assert isinstance(other, MonomialExpansion)

        res = MonomialExpansion({})

        for order, exp in self.expansion.items():
            if order == 0:
                res.expansion[order] = {}

                for monome, coeff in exp.items():
                    if coeff.real != 0:
                        res.expansion[order][monome] = coeff

                        break

                continue

            res.expansion[order] = exp.copy()

            if other.expansion.get(order) is not None:
                for monome, coeff in other.expansion[order].items():
                    if res.expansion[order].get(monome) is None:
                        res.expansion[order][monome] = coeff
                    else:
                        res.expansion[order][monome] += coeff

        for order in other.expansion:
            if res.expansion.get(order) is None:
                res.expansion[order] = other.expansion[order].copy()

        return res

    def extract_order(self, order: int):
        assert order >= 0

        if self.expansion.get(order) is None:
            return MonomialExpansion({})

        return MonomialExpansion({order: self.expansion[order].copy()})

    def reduce(self, monome: Monome):
        newexp = MonomialExpansion({})

        for order, exp in self.expansion.items():
            for m, coeff in exp.items():
                totalorder = order * m.weight()

                if order != 0 and (totalorder % monome.weight() != 0 or totalorder < monome.weight()):
                    continue

                neworder = totalorder // monome.weight()

                if newexp.expansion.get(neworder) is None:
                    newexp.expansion[neworder] = {}

                newexp.expansion[neworder][monome] = coeff

        return newexp

    def up_to_order(self, max_order: int):
        return MonomialExpansion({order: exp for order, exp in self.expansion.items() if order <= max_order})

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
                if m.weight() < 0:
                    m.complex_conjugate = True
                else:
                    m.complex_conjugate = False


        if remove_cc and m.complex_conjugate == True:
            m = m.conjugate()

        if monoms[m] == 0:
            monoms[m] = 1

    return list(monoms.elements())

# def find_fundamental_invariants(variables: list[Variable], n: int, min_order: int = 1, max_order: int = None, remove_cc: bool = True) -> tuple[list[ComplexInvariant], list[Monome]]:
#     if max_order is None:
#         max_order = n

#     fundamentals = []

#     for order in range(min_order, max_order + 1):
#         monoms = generate_monoms(variables, order, n, remove_factorizable=False, remove_cc=remove_cc)

#         for m in monoms:
#             if not m.is_factorisable(n):
#                 continue

#             if not try_to_factorize(m, fundamentals):
#                 fundamentals.append(m)

#     return ([ComplexInvariant(finv, finv.is_real()) for finv in fundamentals], [Monome(finv.variables, complex_conjugate=False, real=False, imag=True) for finv in fundamentals if not finv.is_real()])


def generate_invariants_and_monoms(variables: list[Variable], n: int, min_order: int = 0, max_order: int = None, remove_cc: bool = True) -> tuple[list[Monome], list[Monome], list[Monome]]:
    """
    Generate all invariants, and returns the additional monoms that can appear alongside the appearing monomials
    """
    if max_order is None:
        max_order = n

    invs = []
    rhos = []
    amonoms = []

    for order in range(1, max_order + 1):
        monoms = generate_monoms(variables, order, n, remove_factorizable=False, remove_cc=remove_cc)

        for m in monoms:
            if try_to_factorize(m, invs) or try_to_factorize(m.conjugate(), invs) or m.weight() > n:
                continue

            if m.is_Cn_invariant(n):
                if m.is_sigman_invariant(n):
                    invs.append(Monome(m.variables, complex_conjugate=m.complex_conjugate, invariant_type=InvariantType.full_invariant())) # add monome as invariant
                else:
                    if not m.is_pure_imag():
                        invs.append(Monome(m.variables, complex_conjugate=m.complex_conjugate, invariant_type=InvariantType.real_invariant())) # real part is always invariant

                    invs.append(Monome(m.variables, complex_conjugate=m.complex_conjugate, invariant_type=InvariantType.pseudo_invariant())) # square of imaginary part is invariant
                    rhos.append(Monome(m.variables, complex_conjugate=m.complex_conjugate, invariant_type=InvariantType(invariant=False, real=False, imag=True))) # imaginary part can now appear in monomial expansion
            else:
                amonoms.append(m)

    return (invs, rhos, amonoms)


# def compute_invariants_and_monoms(variables: list[Variable], n: int, min_order: int = 1, max_order: int = None, remove_cc: bool = True) -> tuple[list[Monome], list[Monome]]:
#     if max_order is None:
#         max_order = n

#     monoms = []
#     invariants = []

#     for order in range(1, max_order + 1):
#         raw_monoms = generate_monoms(variables, order, n, remove_factorizable=False, remove_cc=remove_cc)

#         for m in raw_monoms:
#             if m.is_Cn_invariant(n): # weight sum is 0
#                 #if not m.is_pure_imag():
#                 invariants.append(Monome(m.variables, m.complex_conjugate, InvariantType.real_invariant()))

#                 if not m.is_sigman_invariant(n): # anti symmetric under sigma_n
#                     invariants.append(Monome(m.variables + m.variables, m.complex_conjugate, InvariantType.real_invariant()))
#                     monoms.append(Monome(m.variables, m.complex_conjugate, InvariantType.pseudo_invariant()))
#             else:
#                 monoms.append(m)

#     return (
#         #[m for m in monoms if not try_to_factorize(m, invariants) and not try_to_factorize(m.conjugate(), invariants) and abs(m.weight()) <= max_order and (not remove_cc or m.weight() >= 0)],
#         monoms,
#         #[inv for inv in invariants if not try_to_factorize(inv, invariants) and not try_to_factorize(inv.conjugate(), invariants) and abs(m.weight()) <= max_order and (not remove_cc or m.weight() >= 0)]
#         invariants
#     )



# def generate_appearing_monoms(variables: list[Variable], n: int, min_order: int = 0, max_order: int = None, remove_cc: bool = True) -> list[Monome]:
#     if max_order is None:
#         max_order = n

#     #avariables = filter_appearing_variables(variables)
#     finvs, finvs_monoms = find_fundamental_invariants(variables, n, min_order=min_order, max_order=max_order, remove_cc=remove_cc)
#     finvs = list(map(lambda x:x.monome, finvs))
#     monoms = []
#     amonoms = []

#     for order in range(min_order, max_order + 1):
#         monoms.extend(generate_monoms(variables, order, n, remove_factorizable=True, remove_cc=remove_cc))

#     for m in monoms:
#         # need to check for max order since a variable doesn't necessarily have a weight of +/- 1
#         if not try_to_factorize(m, finvs) and not try_to_factorize(m.conjugate(), finvs) and abs(m.weight()) <= max_order:
#             if not remove_cc or m.weight() >= 0:
#                 amonoms.append(m)

#     return amonoms

