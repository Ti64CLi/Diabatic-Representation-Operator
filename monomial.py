#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 21 09:51:22 2025

@author: elidumont
"""

import numpy as np

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

    def order(self) -> int:
        return self.sign * self.symmetry

class Rho:
    def __init__(self, variables: [Variable], orders: [int]):
        assert len(variables) == len(orders)

        self.variables = np.array(variables)
        self.orders = np.array(orders)

    def __repr__(self) -> str:
        if np.all(np.array(self.orders) == 0):
            return "1"

        s = "ρ"

        for i, variable in enumerate(self.variables):
            s += sign2sub(variable.sign)
            s += num2sub(self.orders[i])

        return s

    def __eq__(self, other):
        assert isinstance(other, Rho)

        return np.all(self.variables == other.variables) and np.all(self.orders == other.orders)

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
        return self.__repr__() == other.__repr__()

    def as_rho(self) -> Rho:
        return Rho(self.variables, self.orders)

    def number_of_variables(self) -> int:
        return len(self.variables)

    def order(self) -> int:
        o = 0

        for i, variable in enumerate(self.variables):
            o += variable.symmetry * variable.sign * self.orders[i]

        return o


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

def insert_element(l: list, e) -> list:
    for elem in l:
        if e == elem:
            return l

    return l + [e]

def merge_list(l1: list, l2: list) -> list:
    for e in l2:
        l1 = insert_element(l1, e)

    return l1

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