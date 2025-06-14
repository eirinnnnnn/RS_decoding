from poly import *
from gf import *
from encoder import K, R, N


def rs_decode(received: list[int], erasures: list[int]) -> list[int]:
    """
    Reed-Solomon decoding with errors and erasures (LSB-first).
    Args:
        received: received codeword (length N)
        erasures: list of erasure positions (indices in 0..N-1)

    Returns:
        Decoded message (length K, LSB-first)
    """
    sigma0 = build_erasure_locator(erasures)
    R_erased = erase_positions(received, erasures)
    S = compute_syndrome(R_erased)
    S0 = modified_syndrome(S, sigma0)
    sigma1, omega = solve_key_equation(S0, len(erasures))
    sigma = combine_locators(sigma0, sigma1)
    error_pos = find_error_positions(sigma)
    error_mag = evaluate_error_magnitudes(sigma, omega, error_pos)
    corrected = apply_error_correction(received, error_mag)
    return corrected[-K:]

def debug_rs_decode(received: list[int], erasures: list[int]) -> list[int]:
    print("\n🧪 BEGIN RS DECODING DEBUG 🧪")

    print(f"🔹 Received codeword (LSB-first):\n{received}")
    # print(f"🔹 Erasure positions:\n{erasures}")

    # Step 1: Build sigma₀(x)
    sigma0 = build_erasure_locator(erasures)
    # print(f"\n🔸 Erasure locator σ₀(x): {sigma0}")

    # Step 2: Apply erasures to get R'(x)
    R_erased = erase_positions(received, erasures)
    # print(f"🔸 Erased codeword R'(x): {R_erased}")

    # Step 3: Compute syndromes S₁..S_R
    S = compute_syndrome(R_erased)
    print(f"🔸 Syndrome S(x): {S}")

    # Step 4: Modified syndrome S₀(x) = σ₀(x) · S(x) mod x^R
    S0 = modified_syndrome(S, sigma0)
    # print(f"🔸 Modified syndrome S₀(x): {S0}")

    # Step 5: Solve key equation
    sigma1, omega = solve_key_equation(S0, len(erasures))
    print(f"🔸 σ₁(x) from EE: {sigma1}")
    print(f"🔸 ω(x) from EE:  {omega}")

    # Step 6: Combine locators
    sigma = combine_locators(sigma0, sigma1)
    # print(f"🔸 Combined error locator σ(x) = σ₀ · σ₁: {sigma}")

    # Step 7: Chien search
    error_pos = find_error_positions(sigma)
    print(f"🔸 Error positions found: {error_pos}")

    # Step 8: Forney evaluation
    error_mag = evaluate_error_magnitudes(sigma, omega, error_pos)
    print(f"🔸 Error magnitudes: {error_mag}")

    # Step 9: Correction
    corrected = apply_error_correction(received, error_mag)
    print(f"✅ Corrected codeword: {corrected}")

    print("🧪 END RS DECODING DEBUG 🧪\n")
    return corrected[-K:]


def build_erasure_locator(erasures: list[int]) -> list[int]:
    """
    Build erasure locator polynomial:
        sigma0(x) = product_{i in erasures} (1 - alpha^i x)
    """
    sigma0 = [1]  
    for i in erasures:
        term = [gf_sub(0, EXP_TABLE[i]), 1] 
        sigma0 = poly_mul(sigma0, term)
    return sigma0


def erase_positions(received: list[int], erasures: list[int]) -> list[int]:
    """
    Return R'(x): received with erasure locations zeroed out
    """
    erased = received[:]
    for i in erasures:
        erased[i] = 0
    return erased

def compute_syndrome(received: list[int]) -> list[int]:
    """
    Compute syndrome vector S_j = sum R_i * alpha^{ij} for j = 1..R
    """
    return [poly_eval(received, EXP_TABLE[j]) for j in range(1, R + 1)]

def modified_syndrome(syndrome: list[int], sigma0: list[int]) -> list[int]:
    """
    Compute modified syndrome: S0(x) = sigma0(x) * S(x) mod x^r
    """
    S_poly = syndrome[:]
    S0 = poly_mul(sigma0, S_poly)
    S0_trimmed = S0[:R]  # mod x^r, keep only terms < R
    return S0_trimmed

def solve_key_equation(s0: list[int], e0: int) -> tuple[list[int], list[int]]:
    """
    Solve S0(x) = omega(x) / sigma1(x) using extended Euclidean algorithm
    """
    r_poly = [0] * R + [1]  
    mu = (R - e0) // 2
    nu = (R + e0 - 1) // 2
    sigma1, omega = extended_euclidean(r_poly, s0, mu, nu)
    
    # Normalize by dividing by v(0)
    sigma1_0 = sigma1[0]
    if sigma1_0 == 0:
        raise ZeroDivisionError("sigma1(0) = 0, cannot normalize")

    sigma1 = [gf_div(c, sigma1_0) for c in sigma1]
    omega = [gf_div(c, sigma1_0) for c in omega]
    return sigma1, omega


def combine_locators(sigma0: list[int], sigma1: list[int]) -> list[int]:
    return poly_mul(sigma0, sigma1)


def find_error_positions(sigma: list[int]) -> list[int]:
    """
    Chien search: find i such that sigma(alpha^{-i}) == 0
    """
    positions = []
    for i in range(N):
        x_inv = EXP_TABLE[(63 - i) % 63]  
        if poly_eval(sigma, x_inv) == 0:
            positions.append(i)
    return positions


# def poly_derivative(p: list[int]) -> list[int]:
#     return [gf_mul(c, i) for i, c in enumerate(p)][1:]
def poly_derivative(p: list[int]) -> list[int]:
    """
    Formal derivative in characteristic-2 (LSB-first).
    Even-degree terms disappear; an odd-degree coefficient σ_i
    becomes the coefficient of x^{i-1}.
    """
    if len(p) < 2:
        return [0]          # derivative of constant poly

    deriv = [0] * (len(p) - 1)
    for i in range(1, len(p)):
        if i & 1:           # keep odd exponents
            deriv[i - 1] = p[i]

    # optional – trim leading zeros for canonical form
    while len(deriv) > 1 and deriv[-1] == 0:
        deriv.pop()
    return deriv


def evaluate_error_magnitudes(sigma: list[int], omega: list[int], positions: list[int]) -> dict[int, int]:
    """
    Compute error magnitudes using Forney’s formula
    Returns dict: {position: magnitude}
    """
    error_vector = {}
    sigma_prime = poly_derivative(sigma)
    # print(f"sigma:{sigma}, sigma':{sigma_prime}")
    for i in positions:
        x_inv = EXP_TABLE[(63 - i) % 63]
        num = poly_eval(omega, x_inv)
        denom = poly_eval(sigma_prime, x_inv)
        if denom == 0:
            raise ZeroDivisionError(f"σ'({x_inv}) = 0 during Forney eval")
        error_vector[i] = gf_sub(0, gf_div(num, denom))
        # print(f"x_inv (α^-{i}): {x_inv}, ω(α^-{i}): {num},σ'(α^-{i}): {denom}, magnitude: {error_vector[i]}") 
    return error_vector


def apply_error_correction(received: list[int], error_vector: dict[int, int]) -> list[int]:
    """
    Subtract error values at specified positions
    """
    corrected = received[:]
    for i, mag in error_vector.items():
        corrected[i] = gf_sub(corrected[i], mag)
    return corrected