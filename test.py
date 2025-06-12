from poly import poly_add, poly_mul, poly_divmod, poly_scale, poly_shift, poly_deg, poly_trim
from gf import gf_add, gf_mul, gf_div

def test_poly_add():
    f = [1, 2, 3]
    g = [1, 0, 3]
    expected = [0, 2, 0]
    result = poly_add(f, g)
    assert result == expected, f"poly_add failed: got {result}, expected {expected}"
    print("✅ poly_add passed")

def test_poly_mul():
    f = [1, 2]  # 1*x + 2
    g = [3, 0, 1]  # 3*x^2 + 0*x + 1
    # naive reference: (1x + 2)(3x^2 + 1) = 3x^3 + 2*3x^2 + 1x + 2
    expected = [
    gf_add(gf_mul(1, 3), 0),         # x^3 term
    gf_add(gf_mul(1, 0), gf_mul(2, 3)),  # x^2 term
    gf_add(gf_mul(1, 1), gf_mul(2, 0)),  # x term
    gf_mul(2, 1)                     # constant
]

    result = poly_mul(f, g)
    assert result == expected, f"poly_mul failed: got {result}, expected {expected}"
    print("✅ poly_mul passed")

def test_poly_scale():
    f = [1, 2, 3]
    scalar = 5
    expected = [gf_mul(1, 5), gf_mul(2, 5), gf_mul(3, 5)]
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
    f = [0, 0, 1, 2, 0]
    expected = [1, 2, 0]
    result = poly_trim(f)
    assert result == expected, f"poly_trim failed: got {result}, expected {expected}"
    print("✅ poly_trim passed")

def test_poly_deg():
    assert poly_deg([0, 0, 3, 2]) == 1
    assert poly_deg([0, 0, 0]) == -1
    assert poly_deg([5]) == 0
    print("✅ poly_deg passed")

def test_poly_divmod():
    # f = (x^3 + 2x^2 + 3x + 1), g = (x + 1)
    f = [1, 2, 3, 1]
    g = [1, 1]
    q, r = poly_divmod(f, g)
    # Check q(x)*g(x) + r(x) == f(x)
    lhs = poly_add(poly_mul(q, g), r)
    print(poly_mul(q, g))
    print(q, g, r, f)
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
    print("\n🎉 All polynomial function tests passed!")

if __name__ == "__main__":
    run_all_tests()
