#################
### utilities ###
#################

sign_sup = "⁻⁰⁺"
sign_sub = "₋₀₊"
num_sup = "⁰¹²³⁴⁵⁶⁷⁸⁹"
num_sub = "₀₁₂₃₄₅₆₇₈₉"

def num2sub(n: int) -> str:
    if n == 0:
        return num_sub[0]

    n = abs(n)
    s = ""

    while n > 0:
        s = num_sub[n % 10] + s
        n //= 10

    return s

def num2sup(n: int) -> str:
    if n == 0:
        return num_sup[0]

    n = abs(n)
    s = ""

    while n > 0:
        s = num_sup[n % 10] + s
        n //= 10

    return s

def sign2sub(n: int) -> str:
    if n > 0:
        return sign_sub[2]
    if n < 0:
        return sign_sub[0]
    return sign_sub[1]

def sign2sup(n: int) -> str:
    if n > 0:
        return sign_sup[2]
    if n < 0:
        return sign_sup[0]
    return sign_sup[1]