from gf import gf_add, gf_sub, gf_mul, gf_div

def poly_add(f: list[int], g: list[int]) -> list[int]:
    """
    Add two polynomials over GF(2^6).
    Input: f, g in MSB-first order
    Output: (f + g) in MSB-first order
    """
    max_len = max(len(f), len(g))
    f = [0]*(max_len - len(f)) + f
    g = [0]*(max_len - len(g)) + g
    return [gf_add(a, b) for a, b in zip(f, g)]

def poly_mul(f: list[int], g: list[int]) -> list[int]:
    """
    Multiply two polynomials over GF(2^6).
    Output: (f * g) in MSB-first order
    """
    result = [0] * (len(f) + len(g) - 1)
    for i, a in enumerate(f):
        for j, b in enumerate(g):
            result[i + j] = gf_add(result[i + j], gf_mul(a, b))
    return result

def poly_scale(p: list[int], scalar: int) -> list[int]:
    """
    Scale a polynomial by a scalar over GF(2^6).
    """
    return [gf_mul(c, scalar) for c in p]

def poly_shift(p: list[int], n: int) -> list[int]:
    """
    Multiply polynomial by x^n (append n zeros to the right).
    """
    return p + [0]*n

def poly_deg(p: list[int]) -> int:
    """
    Return the degree of the polynomial.
    Leading zeros are ignored.
    """
    for i, c in enumerate(p):
        if c != 0:
            return len(p) - 1 - i
    return -1  # p(x) = 0

def poly_trim(p: list[int]) -> list[int]:
    """
    Remove leading zeros (return canonical form).
    """
    i = 0
    while i < len(p) and p[i] == 0:
        i += 1
    return p[i:] if i < len(p) else [0]

def poly_divmod(f: list[int], g: list[int]) -> tuple[list[int], list[int]]:
    """
    Divide f(x) by g(x) over GF(2^6), all in MSB-first.
    Returns (quotient, remainder).
    """
    f = poly_trim(f)
    g = poly_trim(g)
    deg_f = poly_deg(f)
    deg_g = poly_deg(g)

    if deg_g < 0:
        raise ZeroDivisionError("division by zero")

    quotient = [0] * (deg_f - deg_g + 1)
    remainder = f[:]

    while poly_deg(remainder) >= deg_g:
        shift = poly_deg(remainder) - deg_g 
        lead_coeff = gf_div(remainder[0], g[0])

        # construct term = lead_coeff * x^shift * g(x)
        scaled_g = poly_scale(g, lead_coeff)
        aligned_g = scaled_g + [0] * shift  # right-padding for MSB alignment

        # quotient[shift-from-left]
        quotient[shift] = lead_coeff

        # update remainder
        remainder = poly_trim(poly_add(remainder, aligned_g))


    return (list(reversed(quotient)), remainder)

def extended_euclidean_limited(a: list[int], b: list[int], mu: int, nu: int) -> tuple[list[int], list[int]]:
    """
    Extended Euclidean Algorithm with early stopping.
    Solves a(x) = x^r, b(x) = modified syndrome
    Stops when deg(σ1) <= mu and deg(ω) <= nu

    Returns:
        - σ1(x) = s(x)
        - ω(x) = r(x)
    """
    r0, r1 = a[:], b[:]
    s0, s1 = [1], [0]

    while True:
        if poly_deg(r1) <= nu and poly_deg(s1) <= mu:
            break
        q, r = poly_divmod(r0, r1)
        r0, r1 = r1, r
        s0, s1 = s1, poly_add(s0, poly_mul(q, s1))

    return s1, r1


#############################################################
def test_extended_euclidean_limited():
    # Example:
    # a(x) = x^5
    # b(x) = 12x^4 + 32x^3 + 5x^2 + 2x + 1 (MSB-first)
    a = [1] + [0]*5             # x^5
    b = [12, 32, 5, 2, 1]       # arbitrary syndrome poly

    # Let's say no erasures: e0 = 0, r = 5
    mu = (5 - 0) // 2   # = 2
    nu = (5 + 0 + 1) // 2 - 1  # = 2

    sigma1, omega = extended_euclidean_limited(a, b, mu, nu)

    # Now check if sigma1 * b ≡ omega mod x^5
    product = poly_mul(sigma1, b)
    _, remainder = poly_divmod(product, a)

    # Check correctness
    if poly_trim(remainder) == poly_trim(omega):
        print("✅ Test passed: sigma1*S0 ≡ omega mod x^r")
    else:
        print("❌ Test failed!")
        print(f"sigma1(x): {sigma1}")
        print(f"omega(x):  {omega}")
        print(f"recomputed omega: {remainder}")



if __name__ == "__main__":
    test_extended_euclidean_limited()
