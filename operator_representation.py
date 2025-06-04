from dataclasses import dataclass
from collections import Counter
from symmetry import Symmetry
from variable import Variable
from monome import Monome
from monomial_expansion import MonomialExpansion
import numpy as np
from utils import *

@dataclass
class Operator:
    # array of monomial expansion
    expansion: np.ndarray[MonomialExpansion]

    def __str__(self) -> str:
        # assuming 2x2 matrices (should be easily generalized)
        s11 = str(self.expansion[0, 0])
        s12 = str(self.expansion[0, 1])
        s21 = str(self.expansion[1, 0])
        s22 = str(self.expansion[1, 1])

        maxw = max(len(s11), len(s12), len(s21), len(s22))

        return f"({s11}{" " * (maxw - len(s11))} | {s12}{" " * (maxw - len(s12))})\n({s21}{" " * (maxw - len(s21))} | {s22}{" " * (maxw - len(s22))})"

    def __add__(self, other):
        assert isinstance(other, Operator)
        assert self.expansion.shape == other.expansion.shape

        n, m = self.expansion.shape
        newop = Operator(self.expansion.copy())

        for i in range(n):
            for j in range(m):
                newop[i, j] += other[i, j]

        return newop

    def __add_matrix(self, monome: Monome, order: int, matrix: np.ndarray):
        assert self.expansion.shape == matrix.shape

        n, m = self.expansion.shape

        for i in range(n):
            for j in range(m):
                self.expansion[i, j] += MonomialExpansion({order: {monome: matrix[i, j]}})

    def __apply_mask(self, mask: np.ndarray):
        assert mask.shape == self.expansion.shape

        n, m = self.expansion.shape

        for i in range(n):
            for j in range(m):
                if mask[i, j] == 0:
                    self.expansion[i, j] = MonomialExpansion({})

    def add_X(self, monome: Monome, order: int, sigma: int, sign: int):
        self.__add_matrix(monome, order, sign * np.array([
            [1, sigma * 1j],
            [sigma * 1j, (-1) ** sigma]
        ]))

    def add_Y(self, monome: Monome, order: int, sigma: int, sign: int):
        self.__add_matrix(monome, order, sign * np.array([
            [1j, -sigma],
            [-sigma, ((-1) ** sigma) * 1j]
        ]))

    def add_X_tilde(self, monome: Monome, order: int, sigma: int, sign: int):
        self.__add_matrix(monome, order, sign * np.array([
            [1, sigma * 1j],
            [-sigma * 1j, 1]
        ]))

    def add_Y_tilde(self, monome: Monome, order: int, sigma: int, sign: int):
        self.__add_matrix(monome, order, sign * np.array([
            [1j, -sigma],
            [sigma, 1j]
        ]))

    def extract_order(self, order: int):
        n, m = self.expansion.shape
        newexp = np.full((n, m), MonomialExpansion({}))

        for i in range(n):
            for j in range(m):
                newexp[i, j] = self.expansion[i, j].extract_order(order)

        return Operator(newexp)

    def apply_states_symmetries(self, n: int, s1: Symmetry, s2: Symmetry):
        if (s1.is_B() or s2.is_B()) and n % 2 != 0:
            raise ValueError("n should be even for a B symmetry")

        mask = np.ones((2, 2))

        if not s1.is_E():
            mask[(s1.value() + 1) % 2, :] = 0
        if not s2.is_E():
            mask[:, (s2.value() + 1) % 2] = 0

        self.__apply_mask(mask)


def A_x(n: int, opsymmetry: Symmetry, s1: Symmetry, s2: Symmetry, max_order: int) -> Operator:
    assert opsymmetry.compute_gamma(n) >= 0
    assert s1.compute_gamma(n) >= 0
    assert s2.compute_gamma(n) >= 0

    variables = [Variable("Q", Symmetry("E", gamma=1))]
    monome = Monome(variables)
    Ax = Operator(np.full((2, 2), MonomialExpansion({})))

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

                    if k > max_order:
                        s += 1
                        continue
                    elif k >= 0:
                        if sigma1 * sigma2 > 0:
                            #print("Adding X^"+str(k)+"_"+str(-sigma2))
                            #Ax += OperatorComponent.X(-sigma2, k, monome)
                            Ax.add_X(monome, k, -sigma2, 1)
                        else:
                            #print("Adding X~^"+str(k)+"_"+str(-sigma2))
                            #Ax += OperatorComponent.X_tilde(-sigma2, k, monome)
                            Ax.add_X_tilde(monome, k, -sigma2, 1)

        if s == 8:
            break

        j += 1

    Ax.apply_states_symmetries(n, s1, s2)

    return Ax

def A_y(n: int, opsymmetry: Symmetry, s1: Symmetry, s2: Symmetry, max_order: int) -> Operator:
    assert opsymmetry.compute_gamma(n) >= 0
    assert s1.compute_gamma(n) >= 0
    assert s2.compute_gamma(n) >= 0

    variables = [Variable("Q", Symmetry("E", gamma=1))]
    monome = Monome(variables)
    Ay = Operator(np.full((2, 2), MonomialExpansion({})))

    if opsymmetry.is_A1() or opsymmetry.is_B1():
        return Ay

    gamma = opsymmetry.compute_gamma(n)
    alpha1 = s1.compute_gamma(n)
    alpha2 = s2.compute_gamma(n)

    j = 0

    while True:
        s = 0

        for sg in [-1, 1]:
            for sigma1 in [-1, 1]:
                for sigma2 in [-1, 1]:
                    k = n * j + gamma + sigma1 * alpha1 - sigma2 * alpha2

                    if k > max_order:
                        s += 1
                        continue
                    elif k >= 0:
                        if sigma1 * sigma2 > 0:
                            #Ay += OperatorComponent.X(-sigma2, k, monome)
                            Ay.add_Y(monome, k, -sigma2, sg)
                        else:
                            #Ay += OperatorComponent.X_tilde(-sigma2, k, monome)
                            Ay.add_Y_tilde(monome, k, -sigma2, sg)

        if s == 8:
            break

        j += 1

    Ay.apply_states_symmetries(n, s1, s2)

    return Ay
