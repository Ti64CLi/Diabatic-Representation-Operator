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
    gamma: Optional[int] = None # only for E_gamma and B1/B2 (if n even, B1/B2 === E_n/2)

    def is_A(self) -> bool:
        return self.label[0] == "A" and self.gamma is None

    def is_A1(self) -> bool:
        return self.label == "A1" and self.gamma is None

    def is_A2(self) -> bool:
        return self.label == "A2" and self.gamma is None

    def is_B(self) -> bool:
        return self.label[0] == "B" and self.gamma is not None

    def is_B1(self) -> bool:
        return self.label == "B1" and self.gamma is not None

    def is_B2(self) -> bool:
        return self.label == "B2" and self.gamma is not None

    def is_E(self) -> bool:
        return self.label == "E" and self.gamma is not None

    def weight(self) -> int:
        return self.gamma if self.is_E() or self.is_B() else 0

    def __str__(self) -> str:
        return f"E_{self.gamma}" if self.is_E() else self.label

    def __eq__(self, other):
        assert isinstance(other, Symmetry)

        return self.label == other.label and self.gamma == other.gamma

    def value(self) -> int:
        values = {"A1": 0, "A2": 1, "B1": 2, "B2": 3, "E": 4}
