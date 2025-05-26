from dataclasses import dataclass
from collections import Counter
from itertools import combinations_with_replacement
from variable import Variable
from symmetry import Symmetry
from utils import *

@dataclass
class Monome:
    variables: list[Variable]

    def weight(self) -> int:
        return sum(v.weight() for v in variables)
    
    def weight_mod(self, n: int) -> int:
        return self.weight() % n
    
    def is_factorisable(self, n: int) -> bool:
        return len(self.variables) > 0 and self.weight_mod(n) == 0
    
    def __str__(self) -> str:
        return ''.join(str(v) for v in self.variables)

def generate_monoms(variables: list[Variable], order: int, n: int, remove_factorizable: bool = True) -> list[Monome]:
    """
    Generates monoms of given order and filters factorizable monoms if needed
    """
    monoms = []

    for combo in combinations_with_replacement(variables, order):
        m = Monome(combo)

        if remove_factorizable and m.is_factorisable(n):
            # add treatement for invariant factors ?
            continue
        
        monoms.append(m)
    
    return m

def generate_variables_list(nvarsym: list[int], n: int) -> list[Variable]:
    """
    Generates a list of variables given the number of variables of each symmetry
    """
    if len(nvarsym) == 0:
        return []
    
    variables = []

    # generate A1 variables
    for i in range(nvarsym[0]):
        variables.append(Variable(f"R{i + 1}", Symmetry("A1")))
    
    # generate A2 variables
    if len(nvarsym) > 1:
        for i in range(nvarsym[1]):
            variables.append(Variable(f"ρ{i + 1}", Symmetry("A2")))
    
    # generate B1 variables
    if len(nvarsym) > 2:
        if nvarsym[2] > 0:
            assert n % 2 == 0

        for i in range(nvarsym[2]):
            variables.append(Variable(f"Rb{i + 1}", Symmetry("B1", gamma=n//2)))
    
    # generate B2 variables
    if len(nvarsym) > 3:
        if nvarsym[3] > 0:
            assert n % 2 == 0

        for i in range(nvarsym[3]):
            variables.append(Variable(f"ρb{i + 1}", Symmetry("B2", gamma=n//2)))
    
    # generate E variables
    if len(nvarsym) > 4:
        for gamma3 in range(4, len(nvarsym)):
            for i in range(nvarsym[gamma3]):
                name = f"Q{i + 1}"
                sym = Symmetry("E", gamma=gamma3 - 3)
                variables.append(Variable(name, sym))
                variables.append(Variable(name, sym, complex_conjugate=True)) # add Q+ and Q-
    
    return variables
