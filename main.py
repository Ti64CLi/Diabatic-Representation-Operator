from operator_representation import *

if __name__ == "__main__":
    H = operator(3, Symmetry("A1"), Symmetry("E", gamma=1), Symmetry("E", gamma=1), [0, 0, 0, 0, 2], 2)
    H11 = H[0, 0]

    # H11[0] is A_x
    print(H11[0].compile2(3, Symmetry("A1"), Symmetry("E", gamma=1), Symmetry("E", gamma=1)))
    print("\n", H11[0])
