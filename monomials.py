#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 16 11:10:43 2025

@author: elidumont
"""

from numpy import dot

def next(l, p, F):
    """
    Returns the next number in base p from l
    
    Args :
        - l : base-p list representation
        - p : base (max exponent)
        - F : number of variables (of p-bits)
    """
    
    i = -1
    l[i] += 1
    
    while l[i] > p:
        if i == -F:
            return []
        
        l[i] = 0
        i -= 1
        l[i] += 1
    
    return l

def monomials_E1(n, p, F):
    """
    Returns the list of monomials of order p
    
    Args :
        - n : rotational invariant (C_nv)
        - p : order of wanted monomials
        - F : number of variables
    """
    
    assert p < n
    
    m_list = []
    m = [0] * F
    
    while True:
        if sum(m) == p:
            m_list.append(m[:])
        
        if m[0] == p:
            break
        
        m = next(m, p, F)
    
    return m_list

def weights(lf, F):
    """
    Computes the weights of each variable given their symmetry (E1, E2, ...)
    
    Args :
        - lf : list of number of variables of each symmetry
        - F : total number of variables (F = sum(lf))
    """
    
    w = [0] * F
    idw = 0
    idf = lf[idw]
    curw = 1
    
    for i in range(F):
        if i == idf:
            idw += 1
            idf += lf[idw]
            curw += 1
        
        w[i] = curw
    
    return w

def next_w(l, p, w, F):
    """
    Computes the next F-uple from l in base p where each bit i is weighted with w[i]
    
    Args :
        - l : base-p previous number
        - p : base (max exponent)
        - w : weight of each bit
        - F : total number of variables (bit length)
    """
    
    i = -1
    l[i] += 1
    
    while l[i] > (p // w[i]):
        if i == -F:
            return []
        
        l[i] = 0
        i -= 1
        l[i] += 1
    
    return l

def monomials_E(n, p, lf):
    """
    Computes every monomials of order p given a number of variables of given symmetries
    
    Args :
        - n : rotational invariant (C_nv)
        - p : order of wanted monomials
        - lf : list of number of variables for each symmetry (E1, E2, ...)
    """
    
    assert p < n
    assert len(lf) <= ((n - 1) // 2)
    
    m_list = []
    F = sum(lf)
    m = [0] * F
    w = weights(lf, F)
    
    while True:
        if dot(m, w) == p:
            m_list.append(m[:])
        
        if m[0] == p:
            break
        
        m = next_w(m, p, w, F)
    
    return m_list

"""
print(monomials_E(10, 5, [2, 1]))
"""
