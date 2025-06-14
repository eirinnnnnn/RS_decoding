from gf import gf_add, gf_sub, gf_mul, gf_div

def poly_add(f: list[int], g: list[int]) -> list[int]:
    """
    Add two polynomials over GF(2^6).
    Input: f, g in LSB-first order
    Output: (f + g) in LSB-first order
    """
    max_len = max(len(f), len(g))
    f =  f + [0]*(max_len - len(f)) 
    g =  g + [0]*(max_len - len(g)) 
    return [gf_add(a, b) for a, b in zip(f, g)]

def poly_mul(f: list[int], g: list[int]) -> list[int]:
    """
    Multiply two polynomials over GF(2^6).
    Output: (f * g) in LSB-first order
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
    return [0]*n + p

def poly_deg(p: list[int]) -> int:
    """
    Return the degree of the polynomial.
    """
    # for i, c in enumerate(p):
    #     if c != 0:
    #         return len(p) - 1 - i
    # return -1  # p(x) = 0
    # p = poly_trim(p)
    # return -1 if len(p) == 1 and p[0] == 0 else len(p)-1 
    for i in reversed(range(len(p))):
        if p[i] != 0:
            return i
    return -1


def poly_trim(p: list[int]) -> list[int]:
    """
    Remove leading zeros (return canonical form).
    """
    i = len(p)-1 
    while i+1 and p[i] == 0:
        i -= 1
    return p[:i+1] if i+1 else [0]


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
        lead_coeff = gf_div(remainder[-1], g[-1])

        # construct term = lead_coeff * x^shift * g(x)
        scaled_g = poly_scale(g, lead_coeff)
        aligned_g =[0] * shift + scaled_g   # left-padding for LSB alignment


        # quotient[shift-from-left]
        quotient[shift] = lead_coeff

        # update remainder
        remainder = poly_trim(poly_add(remainder, aligned_g))


    return (quotient, remainder)
    # return (list(reversed(quotient)), remainder)


def extended_euclidean(a: list[int], b: list[int], mu: int, nu: int) -> tuple[list[int], list[int]]:
    """
    Full Extended Euclidean Algorithm (LSB-first) for polynomials over GF(2^6).
    Follows notation from classical coding theory texts (Section 9.4):
        u(x)a(x) + v(x)b(x) = d(x)

    Args:
        a: Dividend polynomial, LSB-first
        b: Divisor polynomial, LSB-first

    Returns:
        (gcd, u(x), v(x)) in LSB-first
    """

    # Initializations: r_{-1}, r_0
    r_prev = poly_trim(a)
    r_curr = poly_trim(b)

    # u_{-1}, u_0
    u_prev = [1]
    u_curr = [0]

    # v_{-1}, v_0
    v_prev = [0]
    v_curr = [1]

    while poly_deg(r_curr) > nu or poly_deg(v_curr) > mu:
        # print(v_curr, r_curr)
        # print(poly_deg(r_curr), poly_deg(v_curr), nu, mu)
        q, r_next = poly_divmod(r_prev, r_curr)

        # Update r
        r_prev, r_curr = r_curr, r_next

        # Update u_i = u_{i-2} + q_i * u_{i-1}
        u_prev, u_curr = u_curr, poly_trim(poly_add(u_prev, poly_mul(q, u_curr)))

        # Update v_i = v_{i-2} + q_i * v_{i-1}
        v_prev, v_curr = v_curr, poly_trim(poly_add(v_prev, poly_mul(q, v_curr)))
        # print(v_curr, r_curr)
        # print(poly_deg(r_curr), poly_deg(v_curr), nu, mu)

    return v_curr, r_curr  