#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 21 09:54:05 2025

@author: elidumont
"""

from dataclasses import dataclass
from typing import Optional

@dataclass(frozen=True)
class Symmetry:
    label: str                  # can be A1, A2, B1, B2 or E
    gamma: Optional[int] = None # only for E_gamma

    def is_E(self) -> bool:
        return self.label == "E" and self.gamma is not None
    
    def weight(self) -> int:
        return self.gamma if self.is_E() else 0
    
    def __str__(self):
        return f"E_{self.gamma}" if self.is_E() else self.label
    
    def __eq__(self, other):
        assert isinstance(other, Symmetry)

        return self.label == other.label and self.gamma == other.gamma
