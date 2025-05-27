from dataclasses import dataclass
from monome import Monome

@dataclass
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
