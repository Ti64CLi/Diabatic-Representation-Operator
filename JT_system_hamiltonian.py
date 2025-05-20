#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 19 14:09:13 2025

@author: elidumont
"""

import monomials as mn
import operators as op

def hamiltonian(n, alpha, p, name = "H"):
    """
    Computes the Hamiltonian for a E_alpha x e_1 JT system to order p
    H is of A1 symmetry so A_y = 0
    """
    
    assert p < n
    assert alpha >= 0 and alpha <= n // 2
    
    H = op.Operator(name, [])
    
    j = 0
    
    while True:
        s = 0
        
        for sigma in [-1, 0, 1]:
            k = n * j - 2 * sigma * alpha
            
            if k < 0:
                continue
            
            if k > p:
                s += 1
                continue
            
            H += op.OperatorComponent.X(sigma, k)
        
        if s == 3:
            break
        
        j += 1
    
    return H

"""print(hamiltonian(6, 1, 4, "H(1, 1)"))
print(hamiltonian(6, 2, 4, "H(2, 2)"))
print(op.A_x(6, 0, 1, 1, 4))
print(op.A_x(6, 0, 2, 2, 4))
print(op.A_x(6, 0, 1, 2, 4))"""
"""print(op.operator("H(1, 1)", op.Symmetry.A1, 6, 1, 1, p=4))
print(op.A_y(5, 1, 1, 2, 4))
print(op.operator("", op.Symmetry.A1, 3, 1, 0, p=2))"""
print(op.operator("H", op.Symmetry.A1, 0, 3, op.Symmetry.E, 1, op.Symmetry.A1, 0, p=4))