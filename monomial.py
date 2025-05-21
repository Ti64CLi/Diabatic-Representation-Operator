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

    def order(self):
        return self.sign * self.symmetry

class Rho:
    def __init__(self, variables: [Variable], orders: [int]):
        assert len(variables) == len(orders)

        self.variables = variables
        self.orders = orders

    def __repr__(self) -> str:
        s = "ρ"

        for i, variable in enumerate(self.variables):
            s += sign2sub(variable.sign)
            s += num2sub(self.orders[i])

        return s

class Monomial:
    def __init__(self, variables: [Variable], orders: [int]):
        assert len(variables) == len(orders)

        self.variables = variables
        self.orders = orders

    def __repr__(self):
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

    def as_rho(self):
        return Rho(self.variables, self.orders)

    def number_of_variables(self):
        return len(self.variables)

    def order(self):
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

def compute_monomials_E(variables: [Variable], order: int, n: int) -> [Monomial]:
    monomials = []
    F = len(variables)
    monom = [0] * F
    w = weights(variables)

    while True:
        if  np.dot(np.array(m) * signs, w) % p == 0:
            m_list.append(m[:])

        if m[0] >= p:
            break

        m = next_w(m, p, w, F)

    return m_list