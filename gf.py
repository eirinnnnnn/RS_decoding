'''
Define the $GF(2^6)$ by the given generator polynomial.
- exp and log table
- defining the binary operations (add, mult, div, inv, exp)
- polynomial evaluator poly_eval()
'''

# gf.py — GF(2^6) arithmetic for RS(63, 42)

# EXP_TABLE[i] = α^i \forall i \in {0,\ldots,62}
EXP_TABLE = [
    1,  2, 4, 8, 16, 
    32, 3, 6, 12, 24, 
    48, 35, 5, 10, 20, 
    40, 19, 38, 15, 30, 
    60, 59, 53, 41, 17, 
    34, 7, 14, 28, 56, 
    51, 37, 9, 18, 36, 
    11, 22, 44, 27, 54, 
    47, 29, 58, 55, 45, 
    25, 50, 39, 13, 26, 
    52, 43, 21, 42, 23, 
    46, 31, 62, 63, 61, 
    57, 49, 33
][:63]

# LOG_TABLE[a] = i  [α^i = a]
LOG_TABLE = [None] * 64
for i, val in enumerate(EXP_TABLE):
    LOG_TABLE[val] = i

# Generator polynomial coefficients for RS(63, 42)
GENERATOR_POLY_COEFFS = [
    58, 62, 59, 7, 35, 58, 63, 47, 51, 6, 33,
    43, 44, 27, 7, 53, 39, 62, 52, 41, 44, 1
]
# MSB
GENERATOR_POLY_COEFFS = list(reversed(GENERATOR_POLY_COEFFS))

# g(alpha^i) values for i = 0..62, directly from project sheet
G_EVAL_LIST = [
    34, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 16, 13, 43, 41, 48, 15, 11, 52, 33, 
    1, 22, 20, 28, 17, 13, 19, 14, 56, 25, 
    45, 10, 46, 45, 8, 12, 14, 44, 14, 17, 
    26, 31, 22, 59, 29, 52, 31, 57, 48, 45, 
    51, 13
]

# Basic field arithmetic
def gf_add(a, b):
    return a ^ b

def gf_sub(a, b):
    return a^b 

def gf_mul(a, b):
    if a == 0 or b == 0:
        return 0
    return EXP_TABLE[(LOG_TABLE[a] + LOG_TABLE[b]) % 63]

def gf_div(a, b):
    if b == 0:
        raise ZeroDivisionError("division by zero in GF(2^6)")
    if a == 0:
        return 0
    return EXP_TABLE[(LOG_TABLE[a] - LOG_TABLE[b]) % 63]

def gf_inv(a):
    if a == 0:
        raise ZeroDivisionError("inverse of 0 in GF(2^6)")
    return EXP_TABLE[(63 - LOG_TABLE[a]) % 63]

def gf_pow(a, n):
    if a == 0:
        return 0
    return EXP_TABLE[(LOG_TABLE[a] * n) % 63]

# Polynomial evaluation at point x in GF(2^6)
def poly_eval(poly_coeffs, x):
    result = 0
    power = 1
    for coeff in reversed(poly_coeffs):  # MSB 
        result = gf_add(gf_mul(coeff, power), result)
        power = gf_mul(power, x)
    return result

# Verification: does g(alpha^i) match expected values?
def verify_generator_polynomial():
    passed = True
    for i in range(63):
        alpha_i = EXP_TABLE[i]
        expected = G_EVAL_LIST[i]
        actual = poly_eval(GENERATOR_POLY_COEFFS, alpha_i)
        if actual != expected:
            print(f"Mismatch at alpha^{i}: expected {expected}, got {actual}")
            passed = False
    if passed:
        print("All g(alpha^i) values match the official table.")


