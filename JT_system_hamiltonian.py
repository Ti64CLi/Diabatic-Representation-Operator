import operators as op
from symmetry import Symmetry

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


# H11 = op.operator_form("H(1,1)", 6, Symmetry("A1"), Symmetry("E", gamma=1), Symmetry("E", gamma=1), p=4)
# print(H11)
# for p in range(5):
#     print("\nAt order p =", p, ":")
#     print(H11[0].extract_order(p))


# H22 = op.operator_form("H(2,2)", 6, Symmetry("A1"), Symmetry("E", gamma=2), Symmetry("E", gamma=2), p=4)
# print(H22)
# print(H22[0].reduce(2))
# for p in range(5):
#     print("\nAt order p =", p, ":")
#     print(H22[0].extract_order(p))

H12 = op.operator_form("H(1,2)", 6, Symmetry("A1"), Symmetry("E", gamma=1), Symmetry("E", gamma=2), p=4)
print(H12)
for p in range(5):
    print("\nAt order p =", p, ":")
    print(H12[0].extract_order(p))

H21 = op.operator_form("H(2, 1)", 6, Symmetry("A1"), Symmetry("E", gamma=2), Symmetry("E", gamma=1), p=4)
print(H21)
for p in range(5):
    print("\nAt order p =", p, ":")
    print(H21[0].extract_order(p))
