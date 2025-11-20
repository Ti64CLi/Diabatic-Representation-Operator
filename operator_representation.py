from dataclasses import dataclass
from collections import Counter
from symmetry import Symmetry
from variable import Variable, generate_variables_list
from monome import Monome
from monomial_expansion import MonomialExpansion, generate_invariants_and_monoms
import numpy as np
from utils import *

@dataclass
class Operator:
    # group order
    group_order: int
    # operator symmetry
    symmetry: Symmetry
    # states symmetries
    states_symmetries: list[Symmetry]
    # array of monomial expansion
    # corresponding to OperatorComponent
    expansion: np.ndarray[MonomialExpansion]

    def __str__(self) -> str:
        # assuming 2x2 matrices (should be easily generalized)
        s11 = str(self.expansion[0, 0])
        s12 = str(self.expansion[0, 1])
        s21 = str(self.expansion[1, 0])
        s22 = str(self.expansion[1, 1])
        
        maxw1 = max(len(s11), len(s21))
        maxw2 = max(len(s12), len(s22))

        return f"({s11}{" " * (maxw1 - len(s11))} | {s12}{" " * (maxw2 - len(s12))})\n({s21}{" " * (maxw1 - len(s21))} | {s22}{" " * (maxw2 - len(s22))})"

    def __add__(self, other):
        assert isinstance(other, Operator)
        assert self.expansion.shape == other.expansion.shape
        assert self.group_order == other.group_order
        assert self.symmetry == other.symmetry
        assert self.states_symmetries == other.states_symmetries

        n, m = self.expansion.shape
        newop = Operator(self.group_order, self.symmetry, self.states_symmetries, self.expansion.copy())

        for i in range(n):
            for j in range(m):
                newop.expansion[i, j] += other.expansion[i, j]

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

        return Operator(self.group_order, self.symmetry, self.states_symmetries, newexp)

    def apply_states_symmetries(self, n: int, s1: Symmetry, s2: Symmetry):
        if (s1.is_B() or s2.is_B()) and n % 2 != 0:
            raise ValueError("n should be even for a B symmetry")

        mask = np.ones((2, 2))

        if not s1.is_E():
            mask[(s1.value() + 1) % 2, :] = 0
        if not s2.is_E():
            mask[:, (s2.value() + 1) % 2] = 0

        self.__apply_mask(mask)

    def up_to_order(self, max_order: int):
        n, m = self.expansion.shape
        newexp = np.full((n, m), MonomialExpansion({}))

        for i in range(n):
            for j in range(m):
                newexp[i, j] = self.expansion[i, j].up_to_order(max_order)

        return Operator(self.group_order, self.symmetry, self.states_symmetries, newexp)

    def reduce(self, monome: Monome):
        n, m = self.expansion.shape
        newexp = np.full((n, m), MonomialExpansion({}))

        for i in range(n):
            for j in range(m):
                newexp[i, j] = self.expansion[i, j].reduce(monome)

        return Operator(self.group_order, self.symmetry, self.states_symmetries, newexp)

    def compile(self) -> str:
        """
        compile Operator data structure into some basic csv format
        the csv is as follow :
        -> CSV file :
        # metada

        ; section 1 : generale info (for redundancy)
        n, op_sym, state1_sym, state2_sym
        ; *_sym are of the Symmetry type which gives the following values :
        ; A1 -> 0, A2 -> 1, B1 -> 2, B2 -> 2 and E_alpha -> 3 + alpha (E_1 would be 4, E_2 5, etc)

        ; section 2 : variables declaration
        n_A1, n_A2, n_B1, n_B2, n_E_alpha, ...

        ; section 3 : pseudo variables
        ; pseudo variables are declared from variables from section 2 using the following format
        : pvar_sym [pvar_idx]
        +/-Re +/-Im
        [var_sym +/-var_idx var_order]+
        ; the + meaning it can repeat more than once (to represent multi variables pseudo variables)

        ; section 4 : components
        [var_sym +/-var_idx var_order]*
        ; the * means it can repeat 0 or more times
        ; then there are the components terms, the first 4 are for the X part of the operator, and the last 4 are for the Y part of the operator
        +/-Re +/-Im, +/-Re +/-Im, +/-Re +/-Im, +/-Re +/-Im,  +/-Re +/-Im, +/-Re +/-Im, +/-Re +/-Im, +/-Re +/-Im
        <---------------------X part---------------------> | <---------------------Y part--------------------->
        ; and repeat for as much components as there are in the expansion
        """

        s = ""

        return s

    def compile2(self, n: int, op_sym: Symmetry, s1: Symmetry, s2: Symmetry, other: 'Operator' = None) -> str:
        """
        Compiles the Operator (and optionally its Y counterpart) into the project's CSV format.
 
        Args:
            n: The group order (e.g. 3 for C3v).
            op_sym: Symmetry of the operator.
            s1: Symmetry of state 1.
            s2: Symmetry of state 2.
            other: The 'Y' component Operator. If None, Y parts are zeroed.
        """

        # --- Helper Functions ---
        def format_complex(c: complex) -> str:
            re_s = f"{c.real:+.8f}" if c.real != 0 else "+0"
            im_s = f"{c.imag:+.8f}" if c.imag != 0 else "+0"
            return f"{re_s} {im_s}"

        def get_var_sym_code(sym: Symmetry) -> int:
            # A1 -> 0, A2 -> 1, B1 -> 2, B2 -> 3, E_gamma -> 3 + gamma
            val = sym.value() 
            if sym.is_E():
                return 3 + sym.gamma
            return val

        def get_var_idx(var: Variable) -> int:
            # Assumes name format like "Q1", "R2" -> returns 1, 2
            # Filter numeric part
            num = "".join(filter(str.isdigit, var.name))
            return int(num) if num else 0

        # --- 1. Consolidate Terms from X (self) and Y (other) ---
        # Map: MonomialTerm -> {'X': np.ndarray(2,2), 'Y': np.ndarray(2,2)}
        term_map = {}

        def collect_terms(op_obj, component_key):
            if op_obj is None: return

            # op_obj.expansion is np.ndarray[MonomialExpansion]
            # print("op_obj.expansion :", type(op_obj.expansion))
            rows, cols = op_obj.expansion.shape

            for i in range(rows):
                for j in range(cols):
                    # expansion is MonomialExpansion
                    expansion = op_obj.expansion[i, j]
                    # print("\texpansion :", type(expansion))

                    if not expansion.expansion: continue # Skip empty
                    
                    # expansion.expansion is dict[order, dict[Monome, complex]]
                    # print("\t\texpansion.expansion :", type([*expansion.expansion][0]), ",", type([*expansion.expansion.values()][0]))
                    for order_dict in expansion.expansion.values():
                        # print("\t\t\torder_dict :", type([*order_dict][0]), ",", type([*order_dict.values()][0]))
                        for mterm, coeff in order_dict.items():
                            if mterm not in term_map:
                                # print("\t\t\t\tmterm, coeff :", type(mterm), ',', type(coeff))
                                term_map[mterm] = {
                                    'X': np.zeros((2, 2), dtype=complex),
                                    'Y': np.zeros((2, 2), dtype=complex)
                                }

                            term_map[mterm][component_key][i, j] += coeff

        collect_terms(self, 'X')
        collect_terms(other, 'Y')

        # --- 2. Analyze Variables ---
        # We need counts for n_A1, n_A2, n_B1, n_B2, n_E_1...
        # We iterate all terms to find the max index used for each symmetry.
        max_indices = {} # Map sym_code -> max_idx

        for mterm in term_map.keys():
            # Combine variables from rho and monome
            all_vars = mterm.variables
            for v in all_vars:
                code = get_var_sym_code(v.symmetry)
                idx = get_var_idx(v)
                if idx > max_indices.get(code, 0):
                    max_indices[code] = idx

        # --- 3. Build CSV String ---
        lines = []

        # Metadata
        lines.append("# metadata")
        lines.append("; section 1 : general info (for redundancy)")
        states_symmetries = []
        for state_symmetry in self.states_symmetries:
            states_symmetries.append(str(get_var_sym_code(state_symmetry)))
        lines.append(f"{self.group_order}, {get_var_sym_code(self.symmetry)}, {', '.join(states_symmetries)}")

        # Variables Declaration
        # Order: n_A1(0), n_A2(1), n_B1(2), n_B2(3), n_E_1(4), ...
        # We need to determine how many E_alphas exist. 
        # Assuming we go up to the max key found in max_indices
        max_sym_code = max(max_indices.keys()) if max_indices else 0
        counts = []
        # Standard symmetries 0-3 + Es starting at 4
        limit = max(4, max_sym_code + 1) 

        # Special handling: The format asks for n_E_alpha. 
        # If code 4 is E_1, code 5 is E_2.
        # We just output the list of counts ordered by code index? 
        # "n_A1, n_A2, n_B1, n_B2, n_E_alpha, ..." implies specific order.
        # We will output counts for indices 0, 1, 2, 3, 4, 5 ... 

        for i in range(limit):
            counts.append(str(max_indices.get(i, 0)))

        lines.append("; section 2 : variables declaration")
        lines.append(", ".join(counts))

        # Pseudo Variables
        lines.append("; section 3 : pseudo variables")



        # Components
        lines.append("; section 4 : components")

        for mterm, mat_dict in term_map.items():
            # 4a. Build the term definition string: [var_sym +/-var_idx var_order]*
            term_parts = []

            # Count variables in this term
            all_vars = mterm.variables
            var_counts = Counter()
            # We need to distinguish variables by (symmetry, index, conjugate)
            # Variable equality checks name/sym/conjugate.
            for v in all_vars:
                var_counts[v] += 1

            for v, count in var_counts.items():
                sym_code = get_var_sym_code(v.symmetry)
                idx = get_var_idx(v)
                # If conjugate, index is negative (based on usual conventions for this format)
                final_idx = -idx if v.complex_conjugate else idx
                term_parts.append(f"{sym_code} {final_idx} {count}")

            term_def = ", ".join(term_parts)

            # 4b. Build coefficients string
            # X part (2x2) -> Y part (2x2)
            # Flatten row-major: 00, 01, 10, 11
            coeffs = []

            # X Components
            mx = mat_dict['X']
            for r in range(2):
                for c in range(2):
                    coeffs.append(format_complex(mx[r, c]))

            # Y Components
            my = mat_dict['Y']
            for r in range(2):
                for c in range(2):
                    coeffs.append(format_complex(my[r, c]))

            coeffs_str = ", ".join(coeffs)

            # Combine
            lines.append(f"{term_def}, {coeffs_str}")

        return "\n".join(lines)

def A_x(n: int, opsymmetry: Symmetry, s1: Symmetry, s2: Symmetry, max_order: int) -> Operator:
    assert opsymmetry.compute_gamma(n) >= 0
    assert s1.compute_gamma(n) >= 0
    assert s2.compute_gamma(n) >= 0

    variables = [Variable("Q", Symmetry("E", gamma=1))]
    monome = Monome(variables)
    Ax = Operator(n, opsymmetry, [s1, s2], np.full((2, 2), MonomialExpansion({})))

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
                            Ax.add_X(monome, k, -sigma2, 1)
                        else:
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
    Ay = Operator(n, opsymmetry, [s1, s2], np.full((2, 2), MonomialExpansion({})))

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
                            Ay.add_Y(monome, k, -sigma2, sg)
                        else:
                            Ay.add_Y_tilde(monome, k, -sigma2, sg)

        if s == 8:
            break

        j += 1

    Ay.apply_states_symmetries(n, s1, s2)

    return Ay

def operator_form(n: int, opsymmetry: Symmetry, s1: Symmetry, s2: Symmetry, max_order: int) -> tuple[Operator, Operator]:
    return (A_x(n, opsymmetry, s1, s2, max_order), A_y(n, opsymmetry, s1, s2, max_order))

def operator(n: int, opsymmetry: Symmetry, s1: Symmetry, s2: Symmetry, nvarsym: list[int], max_order: int) -> np.ndarray[tuple[Operator, Operator]]:
    """
    Computes the expansion (to order p) of an operator given its symmetry, the symmetry of each state and the symmetry of each variable

    Args :
        - n : type of point group (C_nv)
        - opsymmetry : operator symmetry (A1/2, B1/2, E)
        - s1 : symmetry of the first state
        - s2 : symmetry of the second state
        - nvarsym : list of number of variables of each symmetry
        - max_order : max order of the expansion
    """
    if (opsymmetry.is_B() or s1.is_B() or s2.is_B()) and n % 2 != 0:
        raise ValueError("n should be even for a B symmetry")

    variables = generate_variables_list(nvarsym, n)
    finvs, rhos, monoms = generate_invariants_and_monoms(variables, n)
    n, m = 2, 2
    states = (s1, s2)

    if s1 == s2:
        n, m = 1, 1
        states = (s1,)

    opforms = np.empty((n, m), dtype=object)

    for i in range(n):
        for j in range(m):
            opforms[i, j] = operator_form(n, opsymmetry, states[i], states[j], max_order)

    """opforms = np.array([
        # 11, 12
        [operator_form(n, opsymmetry, s1, s1, max_order), operator_form(n, opsymmetry, s1, s2, max_order)],
        # 21, 22
        [operator_form(n, opsymmetry, s2, s1, max_order), operator_form(n, opsymmetry, s2, s2, max_order)]
    ])"""

    op = np.full((n, m, 2), Operator(n, opsymmetry, [s1, s2], np.full((2, 2), MonomialExpansion({}))))

    for monome in monoms:
        for i in range(n):
            for j in range(m):
                op[i, j][0] += opforms[i, j][0].reduce(monome)
                op[i, j][1] += opforms[i, j][1].reduce(monome)

    return op
