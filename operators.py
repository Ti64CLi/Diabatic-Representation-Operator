#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 19 10:50:49 2025

@author: elidumont
"""

import numpy as np
from symmetry import Symmetry

def sign(n):
    if n > 0:
        return "+"
    elif n < 0:
        return "-"

    return "0"

class ComponentType:
    X = "X"
    Y = "Y"
    X_tilde = "X~"
    Y_tilde = "Y~"

class OperatorComponent:
    """
    Representation of operators components (mainly X^k_sigma and Y^k_sigma) as matrices (2x2 by blocks and function of Q variables)
    """

    def __init__(self, m, k):
        #assert isinstance(ctype, str)
        assert m.shape == (2, 2)
        assert isinstance(k, int)
        assert k >= 0

        #self.ctype = ctype
        self.matrix = m[:, :]
        self.k = k

    def __repr__(self):
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

    def reduce(self, beta):
        assert isinstance(beta, int)
        assert beta > 0
        assert self.k // beta == self.k / beta

        return OperatorComponent(self.matrix, self.k // beta)

    def X(sigma, k):
        assert sigma == 0 or sigma == 1 or sigma == -1
        assert isinstance(k, int)
        assert k >= 0

        return OperatorComponent(np.array([[1, sigma * 1j], [sigma * 1j, (-1)**sigma]]), k)

    def Y(sigma, k):
        assert sigma == 0 or sigma == 1 or sigma == -1
        assert isinstance(k, int)
        assert k >= 0

        return OperatorComponent(np.array([[1j, -sigma], [-sigma, (-1)**sigma]]), k)

    def X_tilde(sigma, k):
        assert sigma != 0 and (sigma == 1 or sigma == -1)
        assert isinstance(k, int)
        assert k >= 0

        return OperatorComponent(np.array([[1, sigma * 1j], [-sigma * 1j, 1]]), k)

    def Y_tilde(sigma, k):
        assert sigma != 0 and (sigma == 1 or sigma == -1)
        assert isinstance(k, int)
        assert k >= 0

        return OperatorComponent(np.array([[1j, -sigma], [sigma, 1j]]), k)

    def apply_symmetry(self, n, s1, alpha1, s2, alpha2):
        if alpha1 == 0 and not s1.is_A():
            raise ValueError("alpha1 should be 0 for a state of A (A1/A2) symmetry")
        if alpha2 == 0 and not s2.is_A():
            raise ValueError("alpha2 should be 0 for a state of A (A1/A2) symmetry")
        if (s1.is_B() or s2.is_B()) and (n // 2 != n / 2):
            raise ValueError("n should be even for a state of B symmetry")
        if (s1.is_B()) and alpha1 != (n // 2):
            raise ValueError("alpha1 should be n/2 for a state of B (B1/B2) symmetry")
        if (s2.is_B()) and alpha2 != (n // 2):
            raise ValueError("alpha2 should be n/2 for a state of B (B1/B2) symmetry")

        if not s1.is_E():
            self.matrix[(s1.value() + 1) % 2, :] = 0
        if not s2.is_E():
            self.matrix[:, (s2.value() + 1) % 2] = 0

class State:
    def __init__(self, symmetry, alpha):
        self.symmetry = symmetry
        self.alpha = alpha

class Operator:
    def __init__(self, name, oc=[]):
        assert isinstance(name, str)

        self.name = name
        self.components = {}

        for component, csign in oc:
            self.__addcomponent(component, csign)

    def __addcomponent(self, component, csign):
        kcomponent = self.components.get(component.k)

        if kcomponent:
            kcomponent, kcsign = kcomponent

            #if np.any(component.matrix != kcomponent.matrix):
            if kcsign == csign:
                kcomponent += component
            else:
                kcomponent -= component

            self.components[component.k] = (kcomponent, kcsign)

            return False

        self.components[component.k] = (component, csign)
        return True

    def __add__(self, other):
        res = Operator(self.name, self.components.values())

        if isinstance(other, OperatorComponent):
            res.__addcomponent(other, 1)

        return res

    def __sub__(self, other):
        res = Operator(self.name, self.components.values())

        if isinstance(other, OperatorComponent):
            res.__addcomponent(other, -1)

        return res

    def __repr__(self):
        s = self.name + " =\n"

        if len(self.components) == 0:
            return s + "0"

        i = 0

        for component, csign in self.components.values():
            if i != 0:
                s += sign(csign)
            #s += " inv*"
            s += str(component)

            i += 1

        return s

    def reduce(self, beta):
        res = Operator(self.name, [])

        for component, csign in self.components.values():
            redc = component.reduce(beta)

            if csign > 0:
                res += redc
            else:
                res -= redc

        return res

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

    def apply_symmetry(self, n, s1, alpha1, s2, alpha2):
        if alpha1 == 0 and not s1.is_A():
            raise ValueError("alpha1 should be 0 for a state of A (A1/A2) symmetry")
        if alpha2 == 0 and not s2.is_A():
            raise ValueError("alpha2 should be 0 for a state of A (A1/A2) symmetry")
        if (s1.is_B() or s2.is_B()) and n % 2 != 0:
            raise ValueError("n should be even for a state of B symmetry")
        if (s1.is_B()) and alpha1 != (n // 2):
            raise ValueError("alpha1 should be n/2 for a state of B (B1/B2) symmetry")
        if (s2.is_B()) and alpha2 != (n // 2):
            raise ValueError("alpha2 should be n/2 for a state of B (B1/B2) symmetry")
        if (s1.is_A()) and alpha1 != 0:
            raise ValueError("alpha1 should be 0 for a state of A (A1/A2) symmetry")
        if (s2.is_A()) and alpha2 != 0:
            raise ValueError("alpha2 should be 0 for a state of A (A1/A2) symmetry")

        mask = np.ones((2, 2))

        if not s1.is_E():
            mask[(s1.value() + 1) % 2, :] = 0
        if not s2.is_E():
            mask[:, (s2.value() + 1) % 2] = 0

        for component, _ in self.components.values():
            component.matrix *= mask

    def extract_order(self, p):
        assert p >= 0

        res = Operator(self.name, [])

        for component, csign in self.components.values():
            if component.k == p:
                res.__addcomponent(component, csign)

        return res

def A_x(n, gamma, alpha1, alpha2, p):
    assert gamma >= 0
    assert alpha1 >= 0
    assert alpha2 >= 0

    Ax = Operator("A_x", [])

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

def A_y(n, gamma, alpha1, alpha2, p):
    assert gamma >= 0
    assert alpha1 >= 0
    assert alpha2 >= 0

    Ay = Operator("A_y", [])

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

def operator(name, opsymmetry, gamma, n, s1, alpha1, s2, alpha2, p=2):
    """
    Computes the expansion (to order p) of an operator given its symmetry and the symmetry of each state

    Args :
        - name : name of the operator (for a display purpose)
        - opsymmetry : operator symmetry (A1/2, B1/2, E)
        - gamma : if the operator is of E symmetry, it's the order of the symmetry (E_gamma), otherwise it has no impact whatsoever
        - n : type of point group (C_nv)
        - s1 : symmetry of the first state
        - alpha1 : the order of the symmetry if the first state is of E symmetry, otherwise is constrained by A1/2 or B1/2 symmetry
        - s2 : symmetry of the second state
        - alpha2 : the order of the symmetry if the second state is of E symmetry, otherwise is constrained by A1/2 or B1/2 symmetry
        - p[=2] : max order of the expansion
    """

    if (opsymmetry.is_B()) and n % 2 != 0:
        raise Exception("n must be even if the operator is of B symmetry")

    op = []

    if opsymmetry.is_A1():
        op = [A_x(n, 0, alpha1, alpha2, p), Operator("A_y")]
    elif opsymmetry.is_A2():
        op = [Operator("A_x"), A_y(n, 0, alpha1, alpha2, p)]
    elif opsymmetry.is_B1():
        op = [A_x(n, n // 2, alpha1, alpha2, p), Operator("A_y")]
    elif opsymmetry.is_B2():
        op = [Operator("A_x"), A_y(n, n // 2, alpha1, alpha2, p)]
    else: # E symmetry
        op = [A_x(n, gamma, alpha1, alpha2, p), A_y(n, gamma, alpha1, alpha2, p)]

    op[0].apply_symmetry(n, s1, alpha1, s2, alpha2)
    op[1].apply_symmetry(n, s1, alpha1, s2, alpha2)

    return op

#def operator(name, opsymmetry, n, )

"""
H = Operator("H(1, 1)", [])

H += OperatorComponent.X(0, 0)
H += OperatorComponent.X(-1, 2)
H += OperatorComponent.X(1, 4)

print(H)
"""