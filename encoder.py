from poly import poly_shift, poly_divmod, poly_add, poly_trim
from gf import GENERATOR_POLY_COEFFS

N = 63  # code length
K = 42  # message length
R = N - K  # number of parity symbols

def rs_encode(message: list[int]) -> list[int]:
    """
    Systematic RS(63, 42) encoder using LSB-first convention.
    Input:
        message: list of K symbols over GF(2^6), LSB-first
    Output:
        codeword: list of N symbols (parity + message), LSB-first
    """
    if len(message) != K:
        raise ValueError(f"Expected message of length {K}, got {len(message)}")

    # Step 1: m(x) * x^r
    m_shifted = poly_shift(message, R)

    # Step 2: Compute remainder: [_, r(x)] = m(x) * x^r mod g(x)
    _, remainder = poly_divmod(m_shifted, GENERATOR_POLY_COEFFS)

    # Step 3: Combine remainder (parity) and message: c(x) = m(x)*x^r - r(x)
    # In LSB-first: [parity..., message...]
    # Pad parity to R symbols (in case deg(r) < R - 1)
    parity = poly_trim(remainder) 

    return (poly_add(parity,  m_shifted))

def test_rs_encode_no_errors():
    msg = [i for i in range(42)]
    codeword = rs_encode(msg)
    assert len(codeword) == 63, "Codeword length mismatch"
    print("✅ RS encode basic test passed")

from encoder import rs_encode
from gf import EXP_TABLE, poly_eval

# Test 1: Check output length
def test_length():
    message = [62-i for i in range(42)]
    codeword = rs_encode(message)
    print(len(codeword))
    assert len(codeword) == 63, "Length of codeword must be 63"
    print("✅ Length test passed.")

# Test 2: Check systematic property (last 42 symbols should equal message)
def test_systematic_property():
    message = [i for i in range(42)]
    codeword = rs_encode(message)
    assert codeword[-42:] == message, "Systematic part (message) incorrect"
    print("✅ Systematic property test passed.")

# Test 3: Check all syndromes of codeword are zero
def test_syndrome_zero():
    message = [i for i in range(42)]
    codeword = rs_encode(message)
    for j in range(1, 22):  # S_1 to S_21
        alpha_j = EXP_TABLE[j % 63]
        val = poly_eval(codeword, alpha_j)
        assert val == 0, f"Syndrome S_{j} = {val} ≠ 0"
    print("✅ All syndromes zero test passed.")



if __name__ == "__main__":
    test_length()
    test_systematic_property()
    test_syndrome_zero()