from encoder import rs_encode, N, K, R
from decoder import debug_rs_decode
from gf import gf_add 
import random

def introduce_errors(codeword, error_positions, error_values):
    corrupted = codeword[:]
    for i, val in zip(error_positions, error_values):
        corrupted[i] = gf_add(corrupted[i], val)
    return corrupted

def introduce_erasures(codeword, erasure_positions):
    corrupted = codeword[:]
    for i in erasure_positions:
        corrupted[i] = 0  # equivalent to ∗ handled by zeroing
    return corrupted

def test_no_errors():
    message = [random.randint(1, 62) for _ in range(K)]
    codeword = rs_encode(message)
    decoded = debug_rs_decode(codeword, [])
    assert decoded == message, "❌ No-error decoding failed"
    print("✅ No-error decoding passed")

def test_with_errors():
    message = [random.randint(1, 62) for _ in range(K)]
    codeword = rs_encode(message)
    error_pos = random.sample(range(N), 5)
    error_vals = [random.randint(1, 62) for _ in error_pos]
    corrupted = introduce_errors(codeword, error_pos, error_vals)
    decoded = debug_rs_decode(corrupted, [])

    print("error position:", error_pos)
    print("error values:", error_vals)

    assert decoded == message, "❌ Error-only decoding failed"
    print("✅ Error-only decoding passed")

def test_with_erasures():
    message = [random.randint(1, 62) for _ in range(K)]
    codeword = rs_encode(message)
    erasure_pos = random.sample(range(N), 5)
    corrupted = introduce_erasures(codeword, erasure_pos)
    decoded = debug_rs_decode(corrupted, erasure_pos)
    assert decoded == message, "❌ Erasure-only decoding failed"
    print("✅ Erasure-only decoding passed")

def test_mixed_errors_erasures():
    message = [random.randint(1, 62) for _ in range(K)]
    codeword = rs_encode(message)
    erasure_pos = random.sample(range(N), 2)
    error_pos = random.sample([i for i in range(N) if i not in erasure_pos], 2)
    error_vals = [random.randint(1, 62) for _ in error_pos]
    corrupted = introduce_errors(codeword, error_pos, error_vals)
    corrupted = introduce_erasures(corrupted, erasure_pos)
    decoded = debug_rs_decode(corrupted, erasure_pos)
    assert decoded == message, "❌ Mixed error-erasure decoding failed"
    print("✅ Mixed error-erasure decoding passed")

if __name__ == "__main__":
    # test_no_errors()
    test_with_errors()
    test_with_erasures()
    test_mixed_errors_erasures()
