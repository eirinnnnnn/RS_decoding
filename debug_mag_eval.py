from typing import List

# Set up the GF(64) primitive field for GF(2^6) using x^6 + x + 1 (0b1000011)
PRIMITIVE_POLY = 0b1000011
GF_SIZE = 64

def generate_gf_tables():
    exp_table = [0] * (2 * GF_SIZE)
    log_table = [0] * GF_SIZE
    x = 1
    for i in range(GF_SIZE - 1):
        exp_table[i] = x
        log_table[x] = i
        x <<= 1
        if x & GF_SIZE:
            x ^= PRIMITIVE_POLY
    for i in range(GF_SIZE - 1, 2 * GF_SIZE):
        exp_table[i] = exp_table[i - (GF_SIZE - 1)]
    return exp_table, log_table

EXP_TABLE, LOG_TABLE = generate_gf_tables()

def gf_add(a: int, b: int) -> int:
    return a ^ b

def gf_mul(a: int, b: int) -> int:
    if a == 0 or b == 0:
        return 0
    return EXP_TABLE[(LOG_TABLE[a] + LOG_TABLE[b]) % (GF_SIZE - 1)]

def gf_div(a: int, b: int) -> int:
    if b == 0:
        raise ZeroDivisionError("Division by zero in GF")
    if a == 0:
        return 0
    return EXP_TABLE[(LOG_TABLE[a] - LOG_TABLE[b]) % (GF_SIZE - 1)]

def poly_eval(p: List[int], x: int) -> int:
    y = 0
    power = 1
    for coeff in p:
        y = gf_add(y, gf_mul(coeff, power))
        power = gf_mul(power, x)
    return y

# Parameters for the case:
omega = [18, 36]           # ω(x) = 37 + 43x
sigma_prime = [44, 39]     # σ'(x) = 19 + 29x
position = 6# Error position index
i = position
# Compute x_inv = α^{-14} = EXP_TABLE[(63 - 14) % 63]
x_inv = EXP_TABLE[(63 - position) % 63]

omega_val = poly_eval(omega, x_inv)
sigma_prime_val = poly_eval(sigma_prime, x_inv)

magnitude = gf_div(omega_val, sigma_prime_val)
print(f"x_inv (α^-{i}): {x_inv}, ω(α^-{i}): {omega_val},σ'(α^-{i}): {sigma_prime_val}, magnitude: {magnitude}")
