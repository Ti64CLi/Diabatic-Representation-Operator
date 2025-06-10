from dataclasses import dataclass
#from monome import Monome

@dataclass
class InvariantType:
    invariant: bool = False
    real: bool = False
    imag: bool = False

    def __str__(self) -> str:
        return str(self.invariant) + str(self.real) + str(self.imag)

    def no_invariant():
        return InvariantType()

    def full_invariant():
        return InvariantType(invariant=True, real=False, imag=False)

    def real_invariant():
        return InvariantType(invariant=True, real=True, imag=False)

    def pseudo_invariant():
        return InvariantType(invariant=True, real=False, imag=True)

    def is_pseudo_invariant(self):
        return self.invariant and self.imag

    def is_real_invariant(self):
        return self.invariant and self.real

    def is_invariant(self):
        return self.invariant


"""@dataclass
class ComplexInvariant:
    monome: Monome
    isreal: bool

    def is_real(self) -> bool:
        return self.isreal

    def imag(self) -> str:
        return f"ρ({str(self.monome)})"

    def real(self) -> str:
        return f"r({str(self.monome)})"

    def __str__(self) -> str:
        if self.is_real():
            return self.real()

        return self.real() + f"+i{self.imag()}"
"""
