#!/usr/bin/env sage
"""
Number field property tests for SageMath comparison.

This script generates deterministic output that can be compared with
the sagemath-ts TypeScript implementation.

Usage:
  sage number_field.sage > ../transcripts/sage/number_field.txt
"""

print("# === Number Field Property Tests ===")
print()

# Set random seed for reproducibility
set_random_seed(12345)

# Test 1: Quadratic field Q(sqrt(2))
print("# === Q(sqrt(2)) ===")
K.<a> = NumberField(x^2 - 2)
print(f"# K.degree()")
print(K.degree())
print()

print(f"# K.discriminant()")
print(K.discriminant())
print()

print(f"# K.signature()")
print(K.signature())
print()

# Element tests
print(f"# a.trace()")
print(a.trace())
print()

print(f"# a.norm()")
print(a.norm())
print()

print(f"# (a+1).trace()")
print((a+1).trace())
print()

print(f"# (a+1).norm()")
print((a+1).norm())
print()

print(f"# a.charpoly()")
print(a.charpoly())
print()

print(f"# a.minpoly()")
print(a.minpoly())
print()

print(f"# (a+1).charpoly()")
print((a+1).charpoly())
print()

# Test 2: Imaginary quadratic field Q(sqrt(-1))
print("# === Q(sqrt(-1)) ===")
K.<i> = NumberField(x^2 + 1)
print(f"# K.degree()")
print(K.degree())
print()

print(f"# K.discriminant()")
print(K.discriminant())
print()

print(f"# K.signature()")
print(K.signature())
print()

print(f"# i.trace()")
print(i.trace())
print()

print(f"# i.norm()")
print(i.norm())
print()

# Test 3: Quadratic field discriminant
print("# === QuadraticField discriminants ===")
for d in [2, 3, 5, -1, -3, -5, -7, -11, 13, 17]:
    K.<a> = QuadraticField(d)
    print(f"# QuadraticField({d}).discriminant()")
    print(K.discriminant())
    print()

# Test 4: Cyclotomic field
print("# === CyclotomicField ===")
for n in [3, 4, 5, 6, 7, 8]:
    K.<z> = CyclotomicField(n)
    print(f"# CyclotomicField({n}).degree()")
    print(K.degree())
    print()

# Test 5: Cubic field
print("# === Cubic field Q(2^(1/3)) ===")
K.<a> = NumberField(x^3 - 2)
print(f"# K.degree()")
print(K.degree())
print()

print(f"# K.signature()")
print(K.signature())
print()

print(f"# a.trace()")
print(a.trace())
print()

print(f"# a.norm()")
print(a.norm())
print()

print(f"# (a^2).trace()")
print((a^2).trace())
print()

print(f"# (a^2).norm()")
print((a^2).norm())
print()

# Test 6: Element arithmetic in Q(sqrt(2))
print("# === Element arithmetic in Q(sqrt(2)) ===")
K.<a> = NumberField(x^2 - 2)

# a^2 = 2
print(f"# a^2")
print(a^2)
print()

# a^3 = 2a
print(f"# a^3")
print(a^3)
print()

# (a+1)^2 = a^2 + 2a + 1 = 3 + 2a
print(f"# (a+1)^2")
print((a+1)^2)
print()

# 1/a = a/2
print(f"# 1/a")
print(1/a)
print()

# is_integral tests
print(f"# a.is_integral()")
print(a.is_integral())
print()

print(f"# (a/2).is_integral()")
print((a/2).is_integral())
print()

print(f"# (a+1).is_integral()")
print((a+1).is_integral())
print()

print("# === End of tests ===")
