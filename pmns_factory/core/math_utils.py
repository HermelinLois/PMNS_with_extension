# ==================================================
# math_utils.py
# Generic math functions for PMNS parameter generation
# ==================================================

def square_and_multiply_no_mod(base, exponent):
    """
    Compute base raised to the power of exponent efficiently.

    Args:
        base : The base to be exponentiated
        exponent : the exponent

    Returns:
        result (same type as base) : value of base ** exponent
    """
    result = 1
    while exponent:
        if exponent & 1:
            result *= base
        base = base ** 2
        exponent >>= 1
    return result


def square_and_multiply_mod(base, exponent, mod):
    """
    Compute base raised to the power of exponent efficiently with mod.

    Args:
        base : The base to be exponentiated
        exponent : the exponent
        mod : the modulus

    Returns:
        result (same type as base) : value of base ** exponent % mod
    """
    resut = 1
    while exponent :
        if exponent&1:
            resut = (resut * base) % mod
        base = pow(base, 2, mod)
        exponent = exponent>>1
    return resut


def square_and_multiply(base, exponent, mod=None):
    """
    Compute base raised to the power of exponent with or without mod.

    Args:
        base : The base to be exponentiated
        exponent : the exponent
        mod : the modulus

    Returns:
        result (same type as base) : value of base ** exponent % mod
    """
    
    if mod is None:
        return square_and_multiply_no_mod(base, exponent)
    return square_and_multiply_mod(base, exponent, mod)

def factorial(n):
    if n <= 1:
        return 1
    return n * factorial(n - 1)