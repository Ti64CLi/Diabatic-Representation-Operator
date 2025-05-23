#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 21 09:51:22 2025

@author: elidumont
"""

import numpy as np

#################
### utilities ###
#################

sign_sup = "⁻⁰⁺"
sign_sub = "₋₀₊"
num_sup = "⁰¹²³⁴⁵⁶⁷⁸⁹"
num_sub = "₀₁₂₃₄₅₆₇₈₉"

def num2sub(n: int) -> str:
    if n == 0:
        return num_sub[0]

    n = abs(n)
    s = ""

    while n > 0:
        s = num_sub[n % 10] + s
        n //= 10

    return s

def num2sup(n: int) -> str:
    if n == 0:
        return num_sup[0]

    n = abs(n)
    s = ""

    while n > 0:
        s = num_sup[n % 10] + s
        n //= 10

    return s

def sign2sub(n: int) -> str:
    return sign_sub[1 + n]

def sign2sup(n: int) -> str:
    return sign_sup[1 + n]

#######################
### Data structures ###
#######################

class Variable:
    def __init__(self, sign, symmetry):
        self.sign = sign
        self.symmetry = symmetry

    def __repr__(self):
        s = "Q" + sign2sub(self.sign)

        if self.symmetry > 1:
            s += num2sub(self.symmetry)

        return s

    def __eq__(self, other):
        assert isinstance(other, Variable)

        return self.sign == other.sign and self.symmetry == other.symmetry

    def __neg__(self):
        return Variable(-self.sign, self.symmetry)

    def order(self) -> int:
        return self.sign * self.symmetry

# could condensate Rho and R class ass one with an extra parameter imaginary for example
# this parameter would determine whether it's a r (i = 0) or rho (i = 1)
class Rho:
    def __init__(self, variables: [Variable], orders: [int], exp=1):
        assert len(variables) == len(orders)

        self.variables = np.array(variables)
        self.orders = np.array(orders)
        self.exp = exp

    def __repr__(self) -> str:
        if np.all(np.array(self.orders) == 0):
            return "1"

        s = "ρ"

        for i, variable in enumerate(self.variables):
            if variable.symmetry != 0:
                s += sign2sub(variable.sign)
                s += num2sub(self.orders[i])
            elif self.orders[i] > 0: # A2 symmetry
                s += "ₐ" + num2sub(i + 1)

        if self.exp > 1:
            s += num2sup(self.exp)

        return s

    def __eq__(self, other):
        assert isinstance(other, Rho)

        return np.all(self.variables == other.variables) and np.all(self.orders == other.orders)

    def as_r(self):
        return R(self.variables, self.orders)

class R:
    def __init__(self, variables: [Variable], orders: [int]):
        assert len(variables) == len(orders)

        self.variables = variables
        self.orders = orders

    def __repr__(self):
        if len(self.variables) == 0:
            return "1"

        s = "r"

        for i, variable in enumerate(self.variables):
            if variable.symmetry != 0:
                s += sign2sub(variable.sign)
                s += num2sub(self.orders[i])
            elif self.orders[i] > 0: # A1 symmetry
                s += "ₐ" + num2sub(i + 1)

        return s

    def as_rho(self):
        return Rho(self.variables, self.orders)

class Monomial:
    def __init__(self, variables: [Variable], orders: [int]):
        assert len(variables) == len(orders)

        self.variables = np.array(variables)
        self.orders = np.array(orders)

    def __repr__(self) -> str:
        if np.all(np.array(self.orders) == 0):
            return "1"

        s = ""

        for i, variable in enumerate(self.variables):
            if self.orders[i] != 0:
                s += "Q" + num2sub(i + 1)
                s += sign2sub(variable.sign)

                if variable.symmetry > 1:
                    s += num2sub(variable.symmetry)

                if self.orders[i] > 1:
                    s += num2sup(self.orders[i])

        return s

    def __eq__(self, other):
        assert isinstance(other, Monomial)
        """
        o1, o2 = [], []

        for i in range(len(self.variables)):
            o1 += [self.variables[i].order() * self.orders[i]]

        for i in range(len(other.variables)):
            o2 += [other.variables[i].order() * other.orders[i]]

        o1.sort(); o2.sort()

        return o1 == o2
        """
        return str(self) == str(other)

    def as_rho(self) -> Rho:
        return Rho(self.variables, self.orders)

    def number_of_variables(self) -> int:
        return len(self.variables)

    def order(self) -> int:
        o = 0

        for i, variable in enumerate(self.variables):
            o += variable.symmetry * variable.sign * self.orders[i]

        return o

class R2(Monomial):
    def __init__(self, variable: int, symmetry: int):
        self.idx = variable
        self.symmetry = symmetry

    def __repr__(self) -> str:
        snum = num2sub(self.idx)
        ssym = ""

        if self.symmetry > 1:
            ssym = num2sub(self.symmetry)

        return "Q" + snum + sign2sub(-1) + ssym + "Q" + snum + sign2sub(1) + ssym

#################################
### General purpose functions ###
#################################

def insert_element(l: list, e) -> list:
    for elem in l:
        if e == elem:
            return l

    return l + [e]

def merge_list(l1: list, l2: list) -> list:
    for e in l2:
        l1 = insert_element(l1, e)

    return l1

def weights(variables: [Variable]) -> [int]:
    w = []

    for variable in variables:
        w += [variable.sign * variable.symmetry]

    return w

def next_orders(orders: [int], w: [int], p: int) -> [int]:
    i = -1
    orders[i] += 1

    while orders[i] > (p // abs(w[i])):
        if i == -len(orders):
            return []

        orders[i] = 0
        i -= 1
        orders[i] += 1

    return orders

def compute_raw_monomials_E(variables: [Variable], order: int, n: int) -> [int]:
    monomials = []
    F = len(variables)
    monom = [0] * F
    w = weights(variables)

    while True:
        relative_order = abs(np.dot(np.array(monom), w))
        if (relative_order == 0 and order % n == 0) or relative_order == order:
            monomials.append(monom[:])

        if monom[0] >= order:
            break

        monom = next_orders(monom, w, order)

    return monomials

def compute_independent_raw_monomials_E(variables: [Variable], order: int, n: int) -> [int]:
    monomials = compute_raw_monomials_E(variables, order, n)
    imonomials = []

    for i, monomi in enumerate(monomials):
        indep = True

        for monomj in monomials[1:i]:
            if all((np.array(monomi) - monomj) >= 0):
                indep = False

                break

        if indep:
            imonomials += [monomi]

    return imonomials

def compute_independent_monomials_E(variables: [Variable], order: int, n: int) -> [Monomial]:
    orders = compute_independent_raw_monomials_E(variables, order, n)
    monomials = []

    for morder in orders:
        monomials += [Monomial(variables, morder)]

    return monomials

###

def compute_independent_rhos_E(variables: [Variable], n: int) -> [Rho]:
    monomials = compute_independent_monomials_E(variables, n, n)

    return [monom.as_rho() for monom in monomials]

def compute_appearing_monomials_E(variables: [Variable], n: int) -> ([Monomial], [Rho]):
    monomials = []
    rhos = compute_independent_rhos_E(variables, n)

    for order in range(n):
        monomials += compute_independent_monomials_E(variables, order, n)

    return (monomials, rhos)

def get_variables(nvarsym: [int], signs: int) -> [Variable]:
    variables = []
    nvar = 0

    for i, nsym in enumerate(nvarsym):
        for j in range(nsym):
            variables += [Variable((-1) ** ((signs >> nvar) & 1), i + 1)]
            nvar += 1

    return variables

def appearing_monomials_E(n: int, nvarsym: [int]) -> [([Monomial], [Rho])]:
    F = sum(nvarsym)
    amonomials = []

    for signs in range(2 ** (F - 1)):
        variables = get_variables(nvarsym, signs)
        #print(variables)
        amonomials += [compute_appearing_monomials_E(variables, n)]

    return amonomials

def merge(amonoms: [([Monomial], [Rho])]) -> ([Monomial], [Rho]):
    if len(amonoms) == 0:
        return [], []

    # add identity beforehand
    monomials = [amonoms[0][0][0]]
    rhos = [amonoms[0][1][0]]

    for amonomials, arhos in amonoms:
        monomials = merge_list(monomials, amonomials)
        #rhos = merge_list(rhos, arhos)
        rhos += arhos[1:] # all rhos should be unique

    return monomials, rhos

def non_trivial_invariants(n: int, amonoms: ([Monomial], [Rho])) -> list:
    monoms, rhos = amonoms

    if len(monoms) == 0 and len(rhos) == 0:
        return []

    inv = []
    variables = []

    if len(monoms) != 0:
        variables = monoms[0].variables
    else:
        variables = rhos[0].variables

    F = len(variables)

    # add Q^n as invariants
    for i, variable in enumerate(variables):
        orders = [0] * F
        orders[i] = n

        inv += [Monomial(variables, orders)]

    # add r and rho^2 as invariants
    for rho in rhos:
        if str(rho) != "1":
            inv += [rho.as_r()]

        inv += [Rho(rho.variables[:], rho.orders[:], 2)]

    return inv

def trivial_invariants(n: int, nvarsym: [int]) -> [R2]:
    tinvariants = []
    nvar = 1

    for sym, nsym in enumerate(nvarsym):
        for i in range(nsym):
            tinvariants += [R2(nvar, sym + 1)]
            nvar += 1

    return tinvariants

def resize_B(n: int, nvarsym: [int]) -> [int]:
    if n % 2 == 0 and len(nvarsym) - 4 < n // 2:
        return nvarsym + [0] * (n // 2 - (len(nvarsym) - 4))

    return nvarsym

def invariants(n: int, nvarsym: [int]) -> list:
    """
    Computes all invariants for a given C_nv symmetry and a list of number of variables of each symmetry with the following order :
        [A1, A2, B1, B2, E, ...]
    It does not yet support B symmetries, but they should regardless be in the nvarsym array to keep track of it
    """

    nvarsym = resize_B(n, nvarsym)

    if n % 2 == 0 and (nvarsym[2] != 0 or nvarsym[3] != 0):
        nvarsym[n // 2] += nvarsym[2] + nvarsym[3]

    nvarEsym = nvarsym[4:]
    invariants = []

    if len(nvarEsym) > 0:
        invariants += trivial_invariants(n, nvarEsym)
        invariants += non_trivial_invariants(n, merge(appearing_monomials_E(n, nvarEsym)))

    # r for A1 symmetry
    a1vars = [Variable(1, 0)] * nvarsym[0]
    for i in range(nvarsym[0]):
        invariants += [R(a1vars, [0] * i + [1] + [0] * (nvarsym[0] - i - 1))]

    # rho² for A2 symmetry
    a2vars = [Variable(1, 0)] * nvarsym[1]

    for i in range(nvarsym[1]):
        invariants += [Rho(a2vars, [0] * i + [1] + [0] * (nvarsym[1] - i - 1), 2)]

    # TODO : Add support for B symmetries (should just be En/2 when n even and none when n odd)

    return invariants





















