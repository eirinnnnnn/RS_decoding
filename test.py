from poly import (
    poly_add, poly_mul, poly_scale, poly_shift,
    poly_deg, poly_trim, poly_divmod, 
)
from gf import gf_add, gf_mul, gf_div


def test_poly_add():
    f = [1, 2, 3]
    g = [1, 0, 3]
    expected = [gf_add(1, 1), gf_add(2, 0), gf_add(3, 3)]
    result = poly_add(f, g)
    assert result == expected, f"poly_add failed: got {result}, expected {expected}"
    print("✅ poly_add passed")

def test_poly_mul():
    f = [1, 2]  # 1 + 2x
    g = [3, 1]  # 3 + x

    expected = [
        gf_mul(1, 3),                                       # x^0
        gf_add(gf_mul(1, 1), gf_mul(2, 3)),                 # x^1
        gf_mul(2, 1)                                        # x^2
    ]

    result = poly_mul(f, g)
    assert result == expected, f"poly_mul failed: got {result}, expected {expected}"
    print("✅ poly_mul passed")

def test_poly_scale():
    f = [1, 2, 3]
    scalar = 5
    expected = [gf_mul(c, scalar) for c in f]
    result = poly_scale(f, scalar)
    assert result == expected, f"poly_scale failed: got {result}, expected {expected}"
    print("✅ poly_scale passed")

def test_poly_shift():
    f = [7, 6]
    expected = [7, 6, 0, 0, 0]
    result = poly_shift(f, 3)
    assert result == expected, f"poly_shift failed: got {result}, expected {expected}"
    print("✅ poly_shift passed")

def test_poly_trim():
    f = [1, 2, 0, 0, 0]
    expected = [1, 2]
    result = poly_trim(f)
    assert result == expected, f"poly_trim failed: got {result}, expected {expected}"
    print("✅ poly_trim passed")

def test_poly_deg():
    assert poly_deg([1, 0, 0]) == 0 
    assert poly_deg([0, 0, 0]) == -1
    assert poly_deg([0, 0, 5]) == 2
    print("✅ poly_deg passed")

def test_poly_divmod():
    f = [1, 2, 3, 1]  # 1 + 2x + 3x^2 + x^3
    g = [1, 1]        # 1 + x
    q, r = poly_divmod(f, g)
    lhs = poly_add(poly_mul(q, g), r)
    assert poly_trim(lhs) == poly_trim(f), f"poly_divmod failed: recomposed {lhs} != {f}"
    print("✅ poly_divmod passed")

def run_all_tests():
    test_poly_add()
    test_poly_mul()
    test_poly_scale()
    test_poly_shift()
    test_poly_trim()
    test_poly_deg()
    test_poly_divmod()
    print("\n🎉 All LSB-first poly tests passed!")

################

from poly import extended_euclidean, poly_trim

def assert_poly_equal(actual, expected, label=""):
    if poly_trim(actual) != poly_trim(expected):
        raise AssertionError(f"❌ {label} failed:\n  Got:      {poly_trim(actual)}\n  Expected: {poly_trim(expected)}")
    print(f"✅ {label} passed")

def test_from_example_9_4():
    # a(x) = x^8
    a = [0]*8 + [1]

    # b(x) = x^6 + x^4 + x^2 + x + 1
    b = [1, 1, 1, 0, 1, 0, 1]

    tests = {
        (0, 7): ([1], [1,1,1,0,1,0,1]),
        (1, 6): ([1], [1,1,1,0,1,0,1]),
        (2, 5): ([1,0,1], [1,1,0,1]),
        (3, 4): ([1,0,1], [1,1,0,1]),
        (4, 3): ([1,0,1], [1,1,0,1]),
        (5, 2): ([0,0,1,1,0,1], [0,0,1]),
        (6, 1): ([1,0,1,1,1,0,1], [1,1]),
        (7, 0): ([1,1,0,1,0,0,1,1], [1])
    }

    for (mu, nu), (expected_sigma, expected_omega) in tests.items():
        sigma, omega = extended_euclidean(a, b, mu, nu)
        label = f"(μ={mu}, ν={nu})"
        assert_poly_equal(sigma, expected_sigma, f"sigma {label}")
        assert_poly_equal(omega, expected_omega, f"omega {label}")

if __name__ == "__main__":
    test_from_example_9_4()
    print("\n🎉 All Euclid tests from Example 9.4 passed.")
