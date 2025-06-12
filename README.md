# RS_decoding



## Structure
```
rs_codec/
├── gf.py       # Galois Field GF(2^6) arithmetic
├── poly.py     # Polynomial arithmetic over GF
├── encoder.py 
├── decoder.py 
├── utils.py   # I/O and testing tools
├── main.py    # Driver script
```

### gf.py: Galois Field Arithmetic for GF(2^6)

This module implements finite field operations over GF(2^6), using the primitive element α defined by the irreducible polynomial:
    $\alpha^6 + \alpha + 1 = 0$

It uses precomputed exponentiation and logarithm tables derived from the project specification.

Features:
- Field operations:
    - gf_add(a, b): Addition (XOR)
    - gf_mul(a, b): Multiplication via log/exp table
    - gf_div(a, b): Division (a / b)
    - gf_inv(a): Multiplicative inverse
    - gf_pow(a, n): Exponentiation
- Polynomial evaluation:
    - poly_eval(poly_coeffs, x): Evaluates a polynomial at x using Horner's method
      (coefficients must be in MSB-first order)
- Generator polynomial verification:
    - verify_generator_polynomial(): Confirms that $g(\alpha^i)$ matches CCC's reference table

Status: Fully tested. All g($\alpha^i$) values match the official table.
