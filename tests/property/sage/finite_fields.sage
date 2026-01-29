"""
Property tests for finite fields - SageMath reference implementation.

This script generates test cases and expected outputs for comparison with sagemath-ts.

Usage:
  sage finite_fields.sage > finite_fields_output.txt
"""

from sage.all import *

# =============================================================================
# GF(p) - Prime Fields Tests
# =============================================================================

print("# ========================================")
print("# Prime Field Tests: GF(p)")
print("# ========================================")

for p in [7, 13, 101, 1009]:
    F = GF(p)
    print(f"# GF({p})")
    print(f"RESULT: order={F.order()}, char={F.characteristic()}")

    # Element operations with fixed elements
    a, b = F(3), F(5)
    print(f"# {a} + {b} in GF({p})")
    print(f"RESULT: {a + b}")
    print(f"# {a} * {b} in GF({p})")
    print(f"RESULT: {a * b}")
    print(f"# {a}^(-1) in GF({p})")
    print(f"RESULT: {a^(-1)}")
    print(f"# {a} - {b} in GF({p})")
    print(f"RESULT: {a - b}")
    print(f"# {a} / {b} in GF({p})")
    print(f"RESULT: {a / b}")
    print(f"# {a}^10 in GF({p})")
    print(f"RESULT: {a^10}")
    print(f"# -{a} in GF({p})")
    print(f"RESULT: {-a}")

    # Multiplicative generator
    g = F.multiplicative_generator()
    print(f"# multiplicative_generator in GF({p})")
    print(f"RESULT: {g}")

    # Generator order
    print(f"# multiplicative_generator order in GF({p})")
    print(f"RESULT: {g.multiplicative_order()}")

# =============================================================================
# GF(p^n) - Extension Fields Tests
# =============================================================================

print("")
print("# ========================================")
print("# Extension Field Tests: GF(p^n)")
print("# ========================================")

extension_fields = [
    (2, 4),   # GF(16)
    (2, 8),   # GF(256)
    (3, 3),   # GF(27)
    (5, 2),   # GF(25)
    (7, 2),   # GF(49)
    (2, 3),   # GF(8)
    (3, 2),   # GF(9)
]

for p, n in extension_fields:
    F = GF(p^n, 'a')
    print(f"# GF({p}^{n})")
    print(f"RESULT: order={F.order()}, char={F.characteristic()}, degree={F.degree()}")

    # Generator and minimal polynomial
    g = F.gen()
    print(f"# gen() in GF({p}^{n})")
    print(f"RESULT: {g}")

    minpoly = g.minimal_polynomial()
    print(f"# gen().minimal_polynomial() in GF({p}^{n})")
    print(f"RESULT: {minpoly}")

    # Conway polynomial (defining polynomial)
    modulus = F.modulus()
    print(f"# modulus() in GF({p}^{n})")
    print(f"RESULT: {modulus}")

    # Generator satisfies modulus
    a = g
    check = modulus(a)
    print(f"# modulus(gen()) should be 0 in GF({p}^{n})")
    print(f"RESULT: {check}")

    # Element arithmetic
    one = F(1)
    print(f"# gen() + 1 in GF({p}^{n})")
    print(f"RESULT: {a + one}")

    print(f"# gen() * gen() in GF({p}^{n})")
    print(f"RESULT: {a * a}")

    print(f"# gen()^{p^n - 1} in GF({p}^{n})")
    print(f"RESULT: {a^(p^n - 1)}")

    # Inverse
    if not a.is_zero():
        print(f"# gen()^(-1) in GF({p}^{n})")
        print(f"RESULT: {a^(-1)}")

        # Verify inverse
        print(f"# gen() * gen()^(-1) in GF({p}^{n})")
        print(f"RESULT: {a * a^(-1)}")

# =============================================================================
# Frobenius Automorphism Tests
# =============================================================================

print("")
print("# ========================================")
print("# Frobenius Automorphism Tests")
print("# ========================================")

for p, n in [(2, 4), (3, 2), (2, 8), (5, 2)]:
    F = GF(p^n, 'a')
    a = F.gen()

    print(f"# Frobenius in GF({p}^{n})")

    # Frobenius: x -> x^p
    frob = a^p
    print(f"# gen()^{p} (Frobenius) in GF({p}^{n})")
    print(f"RESULT: {frob}")

    # Frobenius^n = identity
    frob_n = a^(p^n)
    print(f"# gen()^({p}^{n}) (Frobenius^{n}) should equal gen() in GF({p}^{n})")
    print(f"RESULT: {frob_n == a}")

    # Frobenius is a field homomorphism
    b = a + F(1)
    print(f"# (gen() + 1)^{p} in GF({p}^{n})")
    print(f"RESULT: {b^p}")

    # Verify homomorphism property: (a+b)^p = a^p + b^p in char p
    homomorphism_check = (a + F(1))^p == a^p + F(1)^p
    print(f"# Frobenius homomorphism: (a+1)^p == a^p + 1^p in GF({p}^{n})")
    print(f"RESULT: {homomorphism_check}")

# =============================================================================
# Trace and Norm Tests
# =============================================================================

print("")
print("# ========================================")
print("# Trace and Norm Tests")
print("# ========================================")

for p, n in [(2, 4), (3, 2), (2, 8), (5, 2), (7, 2)]:
    F = GF(p^n, 'a')
    a = F.gen()

    print(f"# Trace and Norm in GF({p}^{n})")

    # Trace: Tr(x) = x + x^p + x^(p^2) + ... + x^(p^(n-1))
    trace_a = a.trace()
    print(f"# trace(gen()) in GF({p}^{n})")
    print(f"RESULT: {trace_a}")

    # Norm: N(x) = x * x^p * x^(p^2) * ... * x^(p^(n-1))
    norm_a = a.norm()
    print(f"# norm(gen()) in GF({p}^{n})")
    print(f"RESULT: {norm_a}")

    # Trace of 1 is n (mod p)
    trace_one = F(1).trace()
    print(f"# trace(1) in GF({p}^{n})")
    print(f"RESULT: {trace_one}")

    # Norm of 1 is 1
    norm_one = F(1).norm()
    print(f"# norm(1) in GF({p}^{n})")
    print(f"RESULT: {norm_one}")

    # Trace is GF(p)-linear: Tr(cx) = c*Tr(x) for c in GF(p)
    c = F.base_ring()(2 % p)
    c_elem = F(c)
    trace_ca = (c_elem * a).trace()
    c_trace_a = c * trace_a
    linear_check = trace_ca == c_trace_a
    print(f"# Trace linearity: Tr(2*gen()) == 2*Tr(gen()) in GF({p}^{n})")
    print(f"RESULT: {linear_check}")

    # Norm is multiplicative: N(xy) = N(x)*N(y)
    b = a + F(1)
    norm_ab = (a * b).norm()
    norm_a_times_norm_b = a.norm() * b.norm()
    mult_check = norm_ab == norm_a_times_norm_b
    print(f"# Norm multiplicativity: N(gen()*(gen()+1)) == N(gen())*N(gen()+1) in GF({p}^{n})")
    print(f"RESULT: {mult_check}")

# =============================================================================
# Conway Polynomial Tests
# =============================================================================

print("")
print("# ========================================")
print("# Conway Polynomial Tests")
print("# ========================================")

# List of fields where Conway polynomials are available
conway_cases = [
    (2, 2), (2, 3), (2, 4), (2, 5), (2, 6), (2, 7), (2, 8),
    (3, 2), (3, 3), (3, 4),
    (5, 2), (5, 3),
    (7, 2), (7, 3),
    (11, 2),
    (13, 2),
]

for p, n in conway_cases:
    try:
        cp = conway_polynomial(p, n)
        print(f"# Conway polynomial C({p},{n})")
        print(f"RESULT: {cp}")

        # Coefficients (excluding leading term)
        coeffs = list(cp)[:-1]  # All coefficients except leading 1
        coeff_str = ",".join(str(c) for c in coeffs)
        print(f"# Conway polynomial C({p},{n}) coefficients")
        print(f"RESULT: [{coeff_str}]")
    except Exception as e:
        print(f"# Conway polynomial C({p},{n})")
        print(f"RESULT: NOT_AVAILABLE")

# =============================================================================
# Element Integer Representation Tests
# =============================================================================

print("")
print("# ========================================")
print("# Element Integer Representation Tests")
print("# ========================================")

for p, n in [(2, 4), (3, 2), (2, 8)]:
    F = GF(p^n, 'a')

    print(f"# Integer representations in GF({p}^{n})")

    # Test a few elements
    for i in range(min(p^n, 10)):
        elem = F.fetch_int(i)
        int_rep = elem.integer_representation()
        print(f"# fetch_int({i}).integer_representation() in GF({p}^{n})")
        print(f"RESULT: {int_rep}")

# =============================================================================
# Field Order Tests (element multiplicative orders)
# =============================================================================

print("")
print("# ========================================")
print("# Multiplicative Order Tests")
print("# ========================================")

for p, n in [(2, 4), (3, 2), (5, 2)]:
    F = GF(p^n, 'a')
    a = F.gen()

    print(f"# Multiplicative orders in GF({p}^{n})")

    # Order of generator
    order_a = a.multiplicative_order()
    print(f"# gen().multiplicative_order() in GF({p}^{n})")
    print(f"RESULT: {order_a}")

    # Order should divide p^n - 1
    divides_check = (p^n - 1) % order_a == 0
    print(f"# multiplicative_order divides {p}^{n}-1 in GF({p}^{n})")
    print(f"RESULT: {divides_check}")

    # Order of 1 is 1
    order_one = F(1).multiplicative_order()
    print(f"# F(1).multiplicative_order() in GF({p}^{n})")
    print(f"RESULT: {order_one}")

# =============================================================================
# Primitive Element Tests
# =============================================================================

print("")
print("# ========================================")
print("# Primitive Element Tests")
print("# ========================================")

for p, n in [(2, 4), (3, 2), (2, 8), (5, 2)]:
    F = GF(p^n, 'a')

    print(f"# Primitive element in GF({p}^{n})")

    # Get primitive element (should have order p^n - 1)
    prim = F.primitive_element()
    prim_order = prim.multiplicative_order()
    expected_order = p^n - 1

    print(f"# primitive_element() in GF({p}^{n})")
    print(f"RESULT: {prim}")

    print(f"# primitive_element().multiplicative_order() in GF({p}^{n})")
    print(f"RESULT: {prim_order}")

    is_primitive = prim_order == expected_order
    print(f"# Is primitive (order == {expected_order})? in GF({p}^{n})")
    print(f"RESULT: {is_primitive}")

# =============================================================================
# Field Embedding Tests (subfield structure)
# =============================================================================

print("")
print("# ========================================")
print("# Subfield Embedding Tests")
print("# ========================================")

# GF(4) embeds in GF(16)
F16 = GF(16, 'a')
F4 = GF(4, 'b')

print("# GF(4) elements in GF(16)")
for i in range(4):
    elem4 = F4.fetch_int(i)
    # Map to GF(16) - elements of GF(4) satisfy x^4 = x
    # They are fixed by Frobenius^2
    print(f"# F4.fetch_int({i}) in GF(4)")
    print(f"RESULT: {elem4}")

# GF(2) embeds in GF(8)
F8 = GF(8, 'a')
F2 = GF(2)

print("# GF(2) elements fixed by Frobenius in GF(8)")
for i in range(2):
    elem = F8(i)
    frob = elem^2
    fixed = frob == elem
    print(f"# F8({i})^2 == F8({i}) in GF(8)")
    print(f"RESULT: {fixed}")

# =============================================================================
# Quadratic Residue Tests
# =============================================================================

print("")
print("# ========================================")
print("# Quadratic Residue Tests (Prime Fields)")
print("# ========================================")

for p in [7, 13, 17, 23]:
    F = GF(p)

    print(f"# Quadratic residues in GF({p})")

    residues = []
    non_residues = []

    for i in range(1, p):
        elem = F(i)
        if elem.is_square():
            residues.append(i)
        else:
            non_residues.append(i)

    print(f"# Quadratic residues in GF({p})")
    print(f"RESULT: {residues}")

    print(f"# Quadratic non-residues in GF({p})")
    print(f"RESULT: {non_residues}")

    # Count should be (p-1)/2 each
    print(f"# Number of residues in GF({p})")
    print(f"RESULT: {len(residues)}")

# =============================================================================
# Square Root Tests
# =============================================================================

print("")
print("# ========================================")
print("# Square Root Tests")
print("# ========================================")

for p in [7, 13, 101]:
    F = GF(p)

    print(f"# Square roots in GF({p})")

    # Test some squares
    for base in [2, 3, 5]:
        elem = F(base)
        sq = elem^2

        if sq.is_square():
            sqrt_val = sq.sqrt()
            # Verify: sqrt^2 == sq
            verify = sqrt_val^2 == sq
            print(f"# sqrt({sq}) in GF({p})")
            print(f"RESULT: {sqrt_val}")
            print(f"# sqrt({sq})^2 == {sq} in GF({p})")
            print(f"RESULT: {verify}")

# =============================================================================
# Field Axioms Verification
# =============================================================================

print("")
print("# ========================================")
print("# Field Axioms Verification")
print("# ========================================")

for p, n in [(2, 4), (3, 2)]:
    F = GF(p^n, 'a')
    a = F.gen()
    b = a + F(1)
    c = a^2

    print(f"# Field axioms in GF({p}^{n})")

    # Commutativity of addition
    print(f"# a + b == b + a")
    print(f"RESULT: {a + b == b + a}")

    # Commutativity of multiplication
    print(f"# a * b == b * a")
    print(f"RESULT: {a * b == b * a}")

    # Associativity of addition
    print(f"# (a + b) + c == a + (b + c)")
    print(f"RESULT: {(a + b) + c == a + (b + c)}")

    # Associativity of multiplication
    print(f"# (a * b) * c == a * (b * c)")
    print(f"RESULT: {(a * b) * c == a * (b * c)}")

    # Distributivity
    print(f"# a * (b + c) == a*b + a*c")
    print(f"RESULT: {a * (b + c) == a*b + a*c}")

    # Additive identity
    print(f"# a + 0 == a")
    print(f"RESULT: {a + F(0) == a}")

    # Multiplicative identity
    print(f"# a * 1 == a")
    print(f"RESULT: {a * F(1) == a}")

    # Additive inverse
    print(f"# a + (-a) == 0")
    print(f"RESULT: {a + (-a) == F(0)}")

    # Multiplicative inverse (for non-zero)
    print(f"# a * a^(-1) == 1")
    print(f"RESULT: {a * a^(-1) == F(1)}")

print("")
print("# ========================================")
print("# All tests completed")
print("# ========================================")
