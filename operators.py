#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 19 10:50:49 2025

@author: elidumont
"""

import numpy as np
from symmetry import Symmetry
from variable import generate_variables_list
from monome import Monome
from monomial_expansion import generate_appearing_monoms, find_fundamental_invariants
from utils import *

class ComponentType:
    X = "X"
    Y = "Y"
    X_tilde = "X~"
    Y_tilde = "Y~"

class OperatorComponent:
    """
    Representation of operators components (mainly X^k_sigma and Y^k_sigma) as matrices (2x2 by blocks and function of Q variables)
    """

    def __init__(self, m: np.ndarray, k: int):
        #assert isinstance(ctype, str)
        assert m.shape == (2, 2)
        assert isinstance(k, int)
        assert k >= 0

        #self.ctype = ctype
        self.matrix = m[:, :]
        self.k = k

    def __str__(self) -> str:
        #s = self.ctype + "^" + str(self.k) + " = \n"
        s = ""

        for i in range(2):
            s += "( "

            for j in range(2):
                if self.matrix[i, j] == 0:
                    s += "     0     "
                elif self.matrix[i, j].imag == 0:
                    s += sign(self.matrix[i, j].real)
                    s += "Re(Q+^" + str(self.k) + ") "
                else:
                    s += sign(self.matrix[i, j].imag)
                    s += "Im(Q+^" + str(self.k) + ") "
            s += ")\n"

        #s += "k = " + str(self.k) + " and sigma = " + str(self.sigma) + "\n"

        return s

    def __add__(self, other):
        assert isinstance(other, OperatorComponent)
        assert self.k == other.k

        return OperatorComponent(self.matrix + other.matrix, self.k)

    def __sub__(self, other):
        assert isinstance(other, OperatorComponent)
        assert self.k == other.k

        return OperatorComponent(self.matrix - other.matrix, self.k)

    def null(k: int):
        return OperatorComponent(np.zeros((2, 2)), k)

    def reduce(self, beta: int):
        assert isinstance(beta, int)
        assert beta > 0

        if self.k % beta != 0:
            return OperatorComponent.null(self.k)

        return OperatorComponent(self.matrix, self.k // beta)

    def X(sigma: int, k: int):
        assert sigma == 0 or sigma == 1 or sigma == -1
        assert isinstance(k, int)
        assert k >= 0

        return OperatorComponent(np.array([[1, sigma * 1j], [sigma * 1j, (-1)**sigma]]), k)

    def Y(sigma: int, k: int):
        assert sigma == 0 or sigma == 1 or sigma == -1
        assert isinstance(k, int)
        assert k >= 0

        return OperatorComponent(np.array([[1j, -sigma], [-sigma, (-1)**sigma]]), k)

    def X_tilde(sigma: int, k: int):
        assert sigma != 0 and (sigma == 1 or sigma == -1)
        assert isinstance(k, int)
        assert k >= 0

        return OperatorComponent(np.array([[1, sigma * 1j], [-sigma * 1j, 1]]), k)

    def Y_tilde(sigma: int, k: int):
        assert sigma != 0 and (sigma == 1 or sigma == -1)
        assert isinstance(k, int)
        assert k >= 0

        return OperatorComponent(np.array([[1j, -sigma], [sigma, 1j]]), k)

    def apply_symmetry(self, n: int, s1: Symmetry, s2: Symmetry):
        if (s1.is_B() or s2.is_B()) and (n // 2 != n / 2):
            raise ValueError("n should be even for a state of B symmetry")

        if not s1.is_E():
            self.matrix[(s1.value() + 1) % 2, :] = 0
        if not s2.is_E():
            self.matrix[:, (s2.value() + 1) % 2] = 0

class Operator:
    def __init__(self, name: str, oc: list[OperatorComponent]=[]):
        assert isinstance(name, str)

        self.name = name
        self.components = {}

        for component, csign in oc:
            self.__addcomponent(component, csign)

    def __addcomponent(self, component: OperatorComponent, csign: int) -> bool:
        kcomponent = self.components.get(component.k)

        if kcomponent is not None:
            kcomponent, kcsign = kcomponent

            if kcsign == csign:
                kcomponent += component
            else:
                kcomponent -= component

            self.components[component.k] = (kcomponent, kcsign)

            return False

        self.components[component.k] = (component, csign)

        return True

    def __add__(self, other):
        assert isinstance(other, OperatorComponent)

        res = Operator(self.name, self.components.values())
        res.__addcomponent(other, 1)

        return res

    def __sub__(self, other):
        assert isinstance(other, OperatorComponent)

        res = Operator(self.name, self.components.values())
        res.__addcomponent(other, -1)

        return res

    def __str__(self) -> str:
        s = self.name + " =\n"

        if len(self.components) == 0:
            return s + "0"

        i = 0

        for component, csign in self.components.values():
            if i != 0:
                s += sign(csign)

            s += str(component)

            i += 1

        return s

    def reduce(self, monome: Monome):
        if monome.weight() == 0:
            return self

        res = Operator(self.name, [])
        beta = monome.weight()

        for component, csign in self.components.values():
            redc = component.reduce(beta)

            if csign > 0:
                res += redc
            else:
                res -= redc

        return res

    # should be redundant
    def explicit(self):
        res = Operator(self.name, [])

        for component, csign in self.components.values():
            if res.components.get(component.k):
                if csign > 0:
                    res.components[component.k] += component
                else:
                    res.components[component.k] -= component
            else:
                res += component

        return res

    def apply_states_symmetries(self, n: int, s1: Symmetry, s2: Symmetry):
        if (s1.is_B() or s2.is_B()) and n % 2 != 0:
            raise ValueError("n should be even for a state of B symmetry")

        mask = np.ones((2, 2))

        if not s1.is_E():
            mask[(s1.value() + 1) % 2, :] = 0
        if not s2.is_E():
            mask[:, (s2.value() + 1) % 2] = 0

        for component, _ in self.components.values():
            component.matrix *= mask

    def extract_order(self, p: int):
        assert p >= 0

        res = Operator(self.name, [])

        for component, csign in self.components.values():
            if component.k == p:
                res.__addcomponent(component, csign)

        return res

def A_x(n: int, opsymmetry: Symmetry, s1: Symmetry, s2: Symmetry, p: int) -> list[Operator]:
    assert opsymmetry.compute_gamma(n) >= 0
    assert s1.compute_gamma(n) >= 0
    assert s2.compute_gamma(n) >= 0

    Ax = Operator("A_x", [])

    if opsymmetry.is_A2() or opsymmetry.is_B2():
        return Ax

    gamma = opsymmetry.compute_gamma(n)
    alpha1 = s1.compute_gamma(n)
    alpha2 = s2.compute_gamma(n)

    j = 0

    while True:
        s = 0

        for sg in [-1, 1]:
            for sigma1 in [-1, 1]:
                for sigma2 in [-1, 1]:
                    k = n * j + sg * gamma + sigma1 * alpha1 + sigma2 * alpha2

                    if k > p:
                        s += 1
                        continue
                    elif k >= 0:
                        if sigma1 * sigma2 > 0:
                            #print("Adding X^"+str(k)+"_"+str(-sigma2))
                            Ax += OperatorComponent.X(-sigma2, k)
                        else:
                            #print("Adding X~^"+str(k)+"_"+str(-sigma2))
                            Ax += OperatorComponent.X_tilde(-sigma2, k)

        if s == 8:
            break

        j += 1

    return Ax

def A_y(n: int, opsymmetry: Symmetry, s1: Symmetry, s2: Symmetry, p: int) -> list[Operator]:
    assert opsymmetry.compute_gamma(n) >= 0
    assert s1.compute_gamma(n) >= 0
    assert s2.compute_gamma(n) >= 0

    Ay = Operator("A_y", [])

    if opsymmetry.is_A1() or opsymmetry.is_B1():
        return Ay

    gamma = opsymmetry.compute_gamma(n)
    alpha1 = s1.compute_gamma(n)
    alpha2 = s2.compute_gamma(n)

    j = 0

    while True:
        s = 0

        for sigma1 in [-1, 1]:
            for sigma2 in [-1, 1]:
                k = n * j + gamma + sigma1 * alpha1 - sigma2 * alpha2

                if k > p:
                    s += 1
                    continue
                elif k >= 0:
                    if sigma1 * sigma2 > 0:
                        Ay += OperatorComponent.X(-sigma2, k)
                    else:
                        Ay += OperatorComponent.X_tilde(-sigma2, k)

        for sigma1 in [-1, 1]:
            for sigma2 in [-1, 1]:
                k = n * j - gamma + sigma1 * alpha1 - sigma2 * alpha2

                if k > p:
                    s += 1
                    continue
                elif k >= 0:
                    if sigma1 * sigma2 > 0:
                        Ay -= OperatorComponent.X(-sigma2, k)
                    else:
                        Ay -= OperatorComponent.X_tilde(-sigma2, k)

        if s == 8:
            break

        j += 1

    return Ay

def operator_form(name: str, n: int, opsymmetry: Symmetry, s1: Symmetry, s2: Symmetry, p=2) -> list[Operator]:
    """
    Computes the expansion (to order p) of an operator given its symmetry and the symmetry of each state

    Args :
        - name : name of the operator (for a display purpose)
        - n : type of point group (C_nv)
        - opsymmetry : operator symmetry (A1/2, B1/2, E)
        - s1 : symmetry of the first state
        - s2 : symmetry of the second state
        - p[=2] : max order of the expansion
    """

    if (opsymmetry.is_B()) and n % 2 != 0:
        raise Exception("n must be even if the operator is of B symmetry")

    op = [A_x(n, opsymmetry, s1, s2, p), A_y(n, opsymmetry, s1, s2, p)]
    op[0].apply_states_symmetries(n, s1, s2)
    op[1].apply_states_symmetries(n, s1, s2)

    return op

def operator(name: str, n: int, opsymmetry: Symmetry, s1: Symmetry, s2: Symmetry, nvarsym: list[int], p=2) -> np.ndarray[list[(Monome, list[Operator])]]:
    variables = generate_variables_list(nvarsym, n)
    monomes = generate_appearing_monoms(variables, n, min_order=1)
    finvs = find_fundamental_invariants(variables, n)

    opform = np.array([
        # 11, 12
        [operator_form(name, n, opsymmetry, s1, s1, p), operator_form(name, n, opsymmetry, s1, s2, p)],
        # 21, 22
        [operator_form(name, n, opsymmetry, s2, s1, p), operator_form(name, n, opsymmetry, s2, s2, p)]
    ])

    op = np.empty((2, 2), dtype=object)

    for i in range(2):
        for j in range(2):
            op[i, j] = []

    for monome in monomes:
        for i in range(2):
            for j in range(2):
                op[i, j].append((monome, [opform[i, j][0].reduce(monome), opform[i, j][1].reduce(monome)]))

    return op

"""
H = Operator("H(1, 1)", [])

H += OperatorComponent.X(0, 0)
H += OperatorComponent.X(-1, 2)
H += OperatorComponent.X(1, 4)

print(H)
"""
