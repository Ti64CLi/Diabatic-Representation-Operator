#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 21 09:51:22 2025

@author: elidumont
"""

from symmetry import Symmetry
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
    if n > 0:
        return sign_sub[2]
    if n < 0:
        return sign_sub[0]
    return sign_sub[1]
    return sign_sub[1 + n]

def sign2sup(n: int) -> str:
    if n > 0:
        return sign_sup[2]
    if n < 0:
        return sign_sup[0]
    return sign_sup[1]
    return sign_sup[1 + n]

#######################
### Data structures ###
#######################

class Variable:
    def __init__(self, symmetry: Symmetry, symmetryorder: int) -> None:
        assert ((symmetry == Symmetry.A1 or symmetry == Symmetry.A2 or symmetry == Symmetry.B1 or symmetry == Symmetry.B2) and symmetryorder == 0) or (symmetry == Symmetry.E and symmetryorder > 0)

        self.symmetry = symmetry
        self.sorder = symmetryorder

    def __repr__(self) -> str:
        s = "Q" + sign2sub(1) + sign2sub(-1)
        s += num2sub(self.sorder)

        return s

    def __eq__(self, other) -> bool:
        assert isinstance(other, Variable)

        return self.symmetry == other.symmetry and self.sorder == other.sorder

    def order(self) -> int:
        return self.sorder

class Monomial:
    """
    Represent the general form of a monomial Q1Q2...Qn

    Members :
        - variables : a list of all possible variables for such a monomial
        - orders : a dict of int * int where the indices are the variables indices and the values are the order. A negative indice indicates a negative variable (Q-)
    """
    def __init__(self, variables: [Variable], orders: {int: int}):
        assert len(orders) <= 2 * len(variables)

        self.variables = np.array(variables)
        self.orders = orders

    def __repr__(self) -> str:
        if len(self.orders) == 0 or all(np.array(list(self.orders.values())) == 0):
            return "1"

        s = ""
        F = len(self.variables)
        indices = [""] * 2 * F

        for i, order in self.orders.items():
            idx = 2 * abs(i) - 1

            if i < 0:
                idx -= 1

            if order > 0:
                indices[idx] = num2sub(abs(i)) + sign2sub(i)

                if self.variables[abs(i) - 1].order() > 1:
                    indices[idx] += num2sub(self.variables[abs(i) - 1].order())

                if order > 1:
                    indices[idx] += num2sup(order)

        for i in indices:
            if len(i) > 0:
                s += "Q" + i

        return s

    def __eq__(self, other):
        assert isinstance(other, Monomial)

        return str(self) == str(other)

    def as_rho(self, order: int):
        return Rho(self.variables, self.orders, order)

    def number_of_variables(self) -> int:
        return len(self.variables)

    def order(self) -> int:
        o = 0

        for i, vorder in self.orders:
            o += self.variables[i].order() * vorder * (i / abs(i))

        return o


# could condensate Rho and R class as one with an extra parameter imaginary for example
# this parameter would determine whether it's a r (i = 0) or rho (i = 1)
class Rho:
    """
    Represent the imaginary part of an invariant QiQj...Qk

    Members :
        - variables : a list of all possible variables Qi
        - orders : an unordered list of tuple representing the order for each variable. The index can be negative when it represent a negative variable Q- (as opposed to Q+). The order must be positive
        - order : order of the rho used. Can be squared for invariant monomial of E symmetry, or 1 for invariant monomial of A2/B2 symmetry
    """
    def __init__(self, variables: [Variable], orders: [(int, int)], order: int):
        assert len(variables) == len(orders)

        self.variables = np.array(variables)
        self.orders = np.array(orders)
        self.order = order

    def __repr__(self) -> str:
        if len(self.orders) == 0 or np.all(self.orders == 0):
            return "1"

        s = "ρ"

        indices = [""] * len(self.variables)

        for i, vorder in self.orders:
            indices[abs(i) - 1] = sign2sub(i) + num2sub(vorder)

        for i in range(len(indices)):
            s += indices[i]

        if self.order > 1:
            s += num2sup(self.order)

        return s

    def __eq__(self, other):
        assert isinstance(other, Rho)

        return str(self) == str(other)

    def as_r(self):
        return R(self.variables, self.orders, self.order)

class R:
    def __init__(self, variables: [Variable], orders: [(int, int)], order: int):
        assert len(variables) == len(orders)

        self.variables = variables
        self.orders = orders
        self.order = order

    def __repr__(self):
        if len(self.variables) == 0:
            return "1"

        s = "r"

        indices = [""] * len(self.variables)

        for i, vorder in self.orders:
            indices[abs(i) - 1] = sign2sub(i) + num2sub(vorder)

        for i in range(len(indices)):
            s += indices[i]

        if self.order > 1:
            s += num2sup(self.order)

        return s

    def as_rho(self):
        return Rho(self.variables, self.orders, self.order)

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

def weights(m: Monomial) -> [int]:
    w = []

    for i, vorder in m.orders.items():
        w += [m.variables[abs(i) - 1].order() * vorder * (i / abs(i))]

    return w

def next_orders(orders: [int], w: [int], pm: int) -> [int]:
    i = -1
    orders[i] += 1

    while orders[i] > (p // abs(w[i])):
        if i == -len(orders):
            return []

        orders[i] = 0
        i -= 1
        orders[i] += 1

    return orders

def compute_raw_monomials_E(m: Monomial, order: int, n: int) -> [int]:
    monomials = []
    F = len(m.variables)
    orders = [0] * F
    w = weights(m)

    while True:
        relative_order = abs(np.dot(np.array(orders), w))
        if (relative_order == 0 and order % n == 0) or relative_order == order:
            monomials.append(orders[:])

        if orders[0] >= order:
            break

        orders = next_orders(orders, w, order)

    return monomials

def compute_independent_raw_monomials_E(m: Monomial, order: int, n: int) -> [int]:
    raw_monomials = compute_raw_monomials_E(m, order, n)
    imonomials = []

    for i, monomi in enumerate(raw_monomials):
        indep = True

        for monomj in raw_monomials[1:i]:
            if all((np.array(monomi) - monomj) >= 0):
                indep = False

                break

        if indep:
            imonomials += [monomi]

    return imonomials

def compute_independent_monomials_E(m: Monomial, order: int, n: int) -> [Monomial]:
    orders = compute_independent_raw_monomials_E(m, order, n)
    monomials = []

    for morder in orders:
        morders = {}

        for i, v in enumerate(morder):
            if m.orders.get(-i):
                morders[-(i + 1)] = v
            else:
                morders[i + 1] = v

        monomials += [Monomial(m.variables, morders)]

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

    # TODO: add Q+Q- as invariants

    # add Q^n as invariants
    for i, variable in enumerate(variables):
        orders = [0] * F
        orders[i] = n

        inv += [Monomial(variables, orders)]

    # add r and rho^2 as invariants
    for rho in rhos:
        if str(rho) != "1":
            inv += [rho.as_r()]

        inv += [Rho(rho.variables[:], rho.orders[:], True)]

    return inv

def trivial_invariants(n: int, nvarsym: [int]) -> [Monomial]:
    tinvariants = []

    for sym, nsym in enumerate(nvarsym):
        for i in range(nsym):
            tinvariants += [Monomial([Variable(1, sym), Variable(-1, sym)], [1, 1])]

    return tinvariants
























