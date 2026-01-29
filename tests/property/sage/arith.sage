# Arithmetic function tests for sagemath-ts comparison
# This file generates test output that TypeScript tests can verify against.
#
# Usage:
#   sage arith.sage > ../transcripts/sage/arith.txt
#
# Or via the comparison runner:
#   ./run-comparison.sh
#
from sage.all import *

# =============================================================================
# GCD Tests
# =============================================================================
print("# gcd(12, 8)")
print(f"RESULT: {gcd(12, 8)}")

print("# gcd(100, 35)")
print(f"RESULT: {gcd(100, 35)}")

print("# gcd(-24, 18)")
print(f"RESULT: {gcd(-24, 18)}")

print("# gcd(0, 5)")
print(f"RESULT: {gcd(0, 5)}")

print("# gcd(1000000, 999999)")
print(f"RESULT: {gcd(1000000, 999999)}")

print("# gcd(2**100, 2**50)")
print(f"RESULT: {gcd(2**100, 2**50)}")

# =============================================================================
# LCM Tests
# =============================================================================
print("# lcm(4, 6)")
print(f"RESULT: {lcm(4, 6)}")

print("# lcm(12, 18)")
print(f"RESULT: {lcm(12, 18)}")

print("# lcm(0, 5)")
print(f"RESULT: {lcm(0, 5)}")

print("# lcm(7, 11)")
print(f"RESULT: {lcm(7, 11)}")

print("# lcm(100, 35)")
print(f"RESULT: {lcm(100, 35)}")

# =============================================================================
# XGCD Tests
# =============================================================================
print("# xgcd(6, 4)")
g, s, t = xgcd(6, 4)
print(f"RESULT: [{g}, {s}, {t}]")

print("# xgcd(240, 46)")
g, s, t = xgcd(240, 46)
print(f"RESULT: [{g}, {s}, {t}]")

print("# xgcd(17, 5)")
g, s, t = xgcd(17, 5)
print(f"RESULT: [{g}, {s}, {t}]")

print("# xgcd(100, 23)")
g, s, t = xgcd(100, 23)
print(f"RESULT: [{g}, {s}, {t}]")

print("# xgcd(-12, 18)")
g, s, t = xgcd(-12, 18)
print(f"RESULT: [{g}, {s}, {t}]")

# =============================================================================
# Factor Tests
# =============================================================================
print("# factor(12)")
f = factor(12)
print(f"RESULT: {list(f)}")

print("# factor(100)")
f = factor(100)
print(f"RESULT: {list(f)}")

print("# factor(1001)")
f = factor(1001)
print(f"RESULT: {list(f)}")

print("# factor(2**31-1)")
f = factor(2**31-1)
print(f"RESULT: {list(f)}")

print("# factor(1)")
f = factor(1)
print(f"RESULT: {list(f)}")

print("# factor(2)")
f = factor(2)
print(f"RESULT: {list(f)}")

print("# factor(1000000)")
f = factor(1000000)
print(f"RESULT: {list(f)}")

print("# factor(97)")
f = factor(97)
print(f"RESULT: {list(f)}")

# =============================================================================
# is_prime Tests
# =============================================================================
print("# is_prime(2)")
print(f"RESULT: {is_prime(2)}")

print("# is_prime(17)")
print(f"RESULT: {is_prime(17)}")

print("# is_prime(100)")
print(f"RESULT: {is_prime(100)}")

print("# is_prime(2**61-1)")
print(f"RESULT: {is_prime(2**61-1)}")

print("# is_prime(1)")
print(f"RESULT: {is_prime(1)}")

print("# is_prime(0)")
print(f"RESULT: {is_prime(0)}")

print("# is_prime(-7)")
print(f"RESULT: {is_prime(-7)}")

print("# is_prime(104729)")
print(f"RESULT: {is_prime(104729)}")

print("# is_prime(561)")
# 561 is a Carmichael number (pseudoprime but not prime)
print(f"RESULT: {is_prime(561)}")

print("# is_prime(2**89-1)")
# Mersenne prime
print(f"RESULT: {is_prime(2**89-1)}")

# =============================================================================
# is_prime_power Tests
# =============================================================================
print("# is_prime_power(8)")
print(f"RESULT: {is_prime_power(8)}")

print("# is_prime_power(9)")
print(f"RESULT: {is_prime_power(9)}")

print("# is_prime_power(12)")
print(f"RESULT: {is_prime_power(12)}")

print("# is_prime_power(17)")
print(f"RESULT: {is_prime_power(17)}")

print("# is_prime_power(1)")
print(f"RESULT: {is_prime_power(1)}")

print("# is_prime_power(8, get_data=True)")
p, k = is_prime_power(8, get_data=True)
print(f"RESULT: [{p}, {k}]")

print("# is_prime_power(81, get_data=True)")
p, k = is_prime_power(81, get_data=True)
print(f"RESULT: [{p}, {k}]")

print("# is_prime_power(12, get_data=True)")
p, k = is_prime_power(12, get_data=True)
print(f"RESULT: [{p}, {k}]")

# =============================================================================
# sqrt_mod Tests (using Mod().sqrt())
# =============================================================================
print("# sqrt_mod(2, 7, all=True)")
r = Mod(2, 7).sqrt(all=True)
print(f"RESULT: {sorted([Integer(x) for x in r]) if r else []}")

print("# sqrt_mod(3, 11, all=True)")
r = Mod(3, 11).sqrt(all=True)
print(f"RESULT: {sorted([Integer(x) for x in r]) if r else []}")

print("# sqrt_mod(4, 13, all=True)")
r = Mod(4, 13).sqrt(all=True)
print(f"RESULT: {sorted([Integer(x) for x in r]) if r else []}")

print("# sqrt_mod(2, 7)")
r = Mod(2, 7).sqrt()
print(f"RESULT: {Integer(r)}")

print("# sqrt_mod(3, 5)")
# 3 is not a QR mod 5
try:
    r = Mod(3, 5).sqrt(extend=False)
    print(f"RESULT: {Integer(r)}")
except ValueError:
    print("RESULT: None")

print("# sqrt_mod(0, 7)")
r = Mod(0, 7).sqrt()
print(f"RESULT: {Integer(r)}")

print("# sqrt_mod(1, 101)")
r = Mod(1, 101).sqrt(all=True)
print(f"RESULT: {sorted([Integer(x) for x in r])}")

print("# sqrt_mod(4, 101)")
r = Mod(4, 101).sqrt(all=True)
print(f"RESULT: {sorted([Integer(x) for x in r])}")

print("# sqrt_mod(2, 17, all=True)")
r = Mod(2, 17).sqrt(all=True)
print(f"RESULT: {sorted([Integer(x) for x in r]) if r else []}")

# =============================================================================
# power_mod Tests
# =============================================================================
print("# power_mod(2, 10, 1000)")
print(f"RESULT: {power_mod(2, 10, 1000)}")

print("# power_mod(3, -1, 7)")
print(f"RESULT: {power_mod(3, -1, 7)}")

print("# power_mod(2, 100, 997)")
print(f"RESULT: {power_mod(2, 100, 997)}")

print("# power_mod(5, 0, 13)")
print(f"RESULT: {power_mod(5, 0, 13)}")

print("# power_mod(7, 1, 101)")
print(f"RESULT: {power_mod(7, 1, 101)}")

print("# power_mod(2, 256, 1000000007)")
print(f"RESULT: {power_mod(2, 256, 1000000007)}")

# =============================================================================
# inverse_mod Tests
# =============================================================================
print("# inverse_mod(3, 7)")
print(f"RESULT: {inverse_mod(3, 7)}")

print("# inverse_mod(5, 11)")
print(f"RESULT: {inverse_mod(5, 11)}")

print("# inverse_mod(17, 101)")
print(f"RESULT: {inverse_mod(17, 101)}")

print("# inverse_mod(2, 1000000007)")
print(f"RESULT: {inverse_mod(2, 1000000007)}")

# =============================================================================
# euler_phi Tests
# =============================================================================
print("# euler_phi(1)")
print(f"RESULT: {euler_phi(1)}")

print("# euler_phi(12)")
print(f"RESULT: {euler_phi(12)}")

print("# euler_phi(100)")
print(f"RESULT: {euler_phi(100)}")

print("# euler_phi(97)")
# Prime: phi(p) = p-1
print(f"RESULT: {euler_phi(97)}")

print("# euler_phi(1000)")
print(f"RESULT: {euler_phi(1000)}")

print("# euler_phi(2**10)")
print(f"RESULT: {euler_phi(2**10)}")

# =============================================================================
# moebius Tests
# =============================================================================
print("# moebius(1)")
print(f"RESULT: {moebius(1)}")

print("# moebius(2)")
print(f"RESULT: {moebius(2)}")

print("# moebius(6)")
# 6 = 2 * 3, two distinct primes
print(f"RESULT: {moebius(6)}")

print("# moebius(4)")
# 4 = 2^2, has squared factor
print(f"RESULT: {moebius(4)}")

print("# moebius(30)")
# 30 = 2 * 3 * 5, three distinct primes
print(f"RESULT: {moebius(30)}")

print("# moebius(5)")
# 5 is prime
print(f"RESULT: {moebius(5)}")

print("# moebius(12)")
# 12 = 4 * 3 = 2^2 * 3, has squared factor
print(f"RESULT: {moebius(12)}")

# =============================================================================
# kronecker_symbol Tests
# =============================================================================
print("# kronecker_symbol(2, 7)")
print(f"RESULT: {kronecker_symbol(2, 7)}")

print("# kronecker_symbol(3, 11)")
print(f"RESULT: {kronecker_symbol(3, 11)}")

print("# kronecker_symbol(5, 13)")
print(f"RESULT: {kronecker_symbol(5, 13)}")

print("# kronecker_symbol(2, 3)")
print(f"RESULT: {kronecker_symbol(2, 3)}")

print("# kronecker_symbol(3, 5)")
print(f"RESULT: {kronecker_symbol(3, 5)}")

print("# kronecker_symbol(0, 5)")
print(f"RESULT: {kronecker_symbol(0, 5)}")

print("# kronecker_symbol(5, 1)")
print(f"RESULT: {kronecker_symbol(5, 1)}")

# =============================================================================
# divisors Tests
# =============================================================================
print("# divisors(12)")
print(f"RESULT: {divisors(12)}")

print("# divisors(100)")
print(f"RESULT: {divisors(100)}")

print("# divisors(1)")
print(f"RESULT: {divisors(1)}")

print("# divisors(97)")
# Prime: only 1 and 97
print(f"RESULT: {divisors(97)}")

print("# divisors(60)")
print(f"RESULT: {divisors(60)}")

# =============================================================================
# sigma Tests (sum of divisors)
# =============================================================================
print("# sigma(12, 1)")
print(f"RESULT: {sigma(12, 1)}")

print("# sigma(12, 0)")
# sigma(n, 0) = number of divisors
print(f"RESULT: {sigma(12, 0)}")

print("# sigma(12, 2)")
# Sum of squares of divisors
print(f"RESULT: {sigma(12, 2)}")

print("# sigma(100, 1)")
print(f"RESULT: {sigma(100, 1)}")

print("# sigma(1, 1)")
print(f"RESULT: {sigma(1, 1)}")

# =============================================================================
# radical Tests
# =============================================================================
print("# radical(12)")
# 12 = 2^2 * 3, radical = 2 * 3 = 6
print(f"RESULT: {radical(12)}")

print("# radical(100)")
# 100 = 2^2 * 5^2, radical = 2 * 5 = 10
print(f"RESULT: {radical(100)}")

print("# radical(1)")
print(f"RESULT: {radical(1)}")

print("# radical(97)")
# Prime: radical = 97
print(f"RESULT: {radical(97)}")

print("# radical(72)")
# 72 = 2^3 * 3^2, radical = 2 * 3 = 6
print(f"RESULT: {radical(72)}")

# =============================================================================
# prime_range Tests
# =============================================================================
print("# prime_range(20)")
print(f"RESULT: {list(prime_range(20))}")

print("# prime_range(10, 30)")
print(f"RESULT: {list(prime_range(10, 30))}")

print("# prime_range(100, 120)")
print(f"RESULT: {list(prime_range(100, 120))}")

print("# prime_range(2, 3)")
print(f"RESULT: {list(prime_range(2, 3))}")

print("# prime_range(1000, 1050)")
print(f"RESULT: {list(prime_range(1000, 1050))}")

# =============================================================================
# valuation Tests
# =============================================================================
print("# valuation(24, 2)")
print(f"RESULT: {valuation(24, 2)}")

print("# valuation(100, 5)")
print(f"RESULT: {valuation(100, 5)}")

print("# valuation(7, 2)")
# 7 is odd, no factors of 2
print(f"RESULT: {valuation(7, 2)}")

print("# valuation(1000, 10)")
print(f"RESULT: {valuation(1000, 10)}")

print("# valuation(81, 3)")
# 81 = 3^4
print(f"RESULT: {valuation(81, 3)}")

print("# valuation(1000000, 2)")
# 1000000 = 2^6 * 5^6
print(f"RESULT: {valuation(1000000, 2)}")

# =============================================================================
# next_prime Tests
# =============================================================================
print("# next_prime(10)")
print(f"RESULT: {next_prime(10)}")

print("# next_prime(97)")
print(f"RESULT: {next_prime(97)}")

print("# next_prime(1)")
print(f"RESULT: {next_prime(1)}")

print("# next_prime(1000)")
print(f"RESULT: {next_prime(1000)}")

print("# next_prime(2**20)")
print(f"RESULT: {next_prime(2**20)}")

# =============================================================================
# previous_prime Tests
# =============================================================================
print("# previous_prime(10)")
print(f"RESULT: {previous_prime(10)}")

print("# previous_prime(100)")
print(f"RESULT: {previous_prime(100)}")

print("# previous_prime(97)")
print(f"RESULT: {previous_prime(97)}")

print("# previous_prime(1000)")
print(f"RESULT: {previous_prime(1000)}")

print("# previous_prime(3)")
print(f"RESULT: {previous_prime(3)}")

# =============================================================================
# isqrt Tests
# =============================================================================
print("# isqrt(100)")
print(f"RESULT: {isqrt(100)}")

print("# isqrt(99)")
print(f"RESULT: {isqrt(99)}")

print("# isqrt(0)")
print(f"RESULT: {isqrt(0)}")

print("# isqrt(1)")
print(f"RESULT: {isqrt(1)}")

print("# isqrt(2**100)")
print(f"RESULT: {isqrt(2**100)}")

# =============================================================================
# is_square Tests
# =============================================================================
print("# is_square(100)")
print(f"RESULT: {is_square(100)}")

print("# is_square(99)")
print(f"RESULT: {is_square(99)}")

print("# is_square(0)")
print(f"RESULT: {is_square(0)}")

print("# is_square(1)")
print(f"RESULT: {is_square(1)}")

print("# is_square(-1)")
print(f"RESULT: {is_square(-1)}")

print("# is_square(2**100)")
print(f"RESULT: {is_square(2**100)}")

# =============================================================================
# is_squarefree Tests
# =============================================================================
print("# is_squarefree(30)")
# 30 = 2 * 3 * 5, squarefree
print(f"RESULT: {is_squarefree(30)}")

print("# is_squarefree(12)")
# 12 = 4 * 3 = 2^2 * 3, not squarefree
print(f"RESULT: {is_squarefree(12)}")

print("# is_squarefree(1)")
print(f"RESULT: {is_squarefree(1)}")

print("# is_squarefree(97)")
# Prime: squarefree
print(f"RESULT: {is_squarefree(97)}")

print("# is_squarefree(100)")
# 100 = 4 * 25 = 2^2 * 5^2, not squarefree
print(f"RESULT: {is_squarefree(100)}")

# =============================================================================
# squarefree_part Tests
# =============================================================================
print("# squarefree_part(12)")
# 12 = 3 * 4 = 3 * 2^2, squarefree part = 3
print(f"RESULT: {squarefree_part(12)}")

print("# squarefree_part(72)")
# 72 = 2 * 36 = 2 * 6^2, squarefree part = 2
print(f"RESULT: {squarefree_part(72)}")

print("# squarefree_part(100)")
# 100 = 1 * 10^2, squarefree part = 1
print(f"RESULT: {squarefree_part(100)}")

print("# squarefree_part(30)")
# 30 = 30 * 1^2, squarefree part = 30
print(f"RESULT: {squarefree_part(30)}")

print("# squarefree_part(1)")
print(f"RESULT: {squarefree_part(1)}")

# =============================================================================
# prime_divisors (prime_factors) Tests
# =============================================================================
print("# prime_divisors(12)")
print(f"RESULT: {prime_divisors(12)}")

print("# prime_divisors(100)")
print(f"RESULT: {prime_divisors(100)}")

print("# prime_divisors(1)")
print(f"RESULT: {prime_divisors(1)}")

print("# prime_divisors(97)")
print(f"RESULT: {prime_divisors(97)}")

print("# prime_divisors(1001)")
print(f"RESULT: {prime_divisors(1001)}")

# =============================================================================
# CRT Tests
# =============================================================================
print("# CRT(2, 3, 3, 5)")
print(f"RESULT: {CRT(2, 3, 3, 5)}")

print("# CRT(1, 2, 3, 7)")
print(f"RESULT: {CRT(1, 2, 3, 7)}")

print("# CRT(0, 0, 5, 7)")
print(f"RESULT: {CRT(0, 0, 5, 7)}")

print("# CRT_list([2, 3, 2], [3, 5, 7])")
print(f"RESULT: {CRT_list([2, 3, 2], [3, 5, 7])}")

# =============================================================================
# number_of_divisors Tests
# =============================================================================
print("# number_of_divisors(12)")
print(f"RESULT: {number_of_divisors(12)}")

print("# number_of_divisors(100)")
print(f"RESULT: {number_of_divisors(100)}")

print("# number_of_divisors(1)")
print(f"RESULT: {number_of_divisors(1)}")

print("# number_of_divisors(97)")
print(f"RESULT: {number_of_divisors(97)}")

# =============================================================================
# Large number tests to stress test implementations
# =============================================================================
print("# gcd(2**64, 2**32)")
print(f"RESULT: {gcd(2**64, 2**32)}")

print("# is_prime(2**127-1)")
# Mersenne prime M127
print(f"RESULT: {is_prime(2**127-1)}")

print("# euler_phi(2**32)")
print(f"RESULT: {euler_phi(2**32)}")

print("# next_prime(2**64)")
print(f"RESULT: {next_prime(2**64)}")

print("")
print("# === END OF TESTS ===")
