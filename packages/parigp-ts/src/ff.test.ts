/**
 * @module ff.test
 * @description Tests for finite field (Fp) operations
 *
 * Test vectors are derived from PARI/GP computations and known mathematical properties.
 */

import { describe, expect, it } from 'bun:test';
import {
  Fp_add,
  Fp_addmul,
  Fp_center,
  Fp_div,
  Fp_double,
  Fp_halve,
  Fp_inv,
  Fp_issquare,
  Fp_mul,
  Fp_neg,
  Fp_pow,
  Fp_red,
  Fp_sqr,
  Fp_sqrt,
  Fp_sub,
  gcd,
  kronecker,
  xgcd,
} from './ff.js';

// Test primes of various sizes
const p7 = 7n; // Small prime
const p11 = 11n; // Small prime
const p23 = 23n; // p = 3 mod 4
const p29 = 29n; // p = 1 mod 4
const p97 = 97n; // p = 1 mod 8
const p1009 = 1009n; // Larger prime, p = 1 mod 8
const p65537 = 65537n; // Fermat prime F_4

describe('Fp_red - Reduction mod p', () => {
  it('should reduce positive integers', () => {
    expect(Fp_red(10n, p7)).toBe(3n);
    expect(Fp_red(7n, p7)).toBe(0n);
    expect(Fp_red(15n, p7)).toBe(1n);
  });

  it('should reduce negative integers', () => {
    expect(Fp_red(-1n, p7)).toBe(6n);
    expect(Fp_red(-7n, p7)).toBe(0n);
    expect(Fp_red(-8n, p7)).toBe(6n);
  });

  it('should handle zero', () => {
    expect(Fp_red(0n, p7)).toBe(0n);
    expect(Fp_red(0n, p1009)).toBe(0n);
  });
});

describe('Fp_add - Addition mod p', () => {
  it('should add without overflow', () => {
    expect(Fp_add(2n, 3n, p7)).toBe(5n);
    expect(Fp_add(0n, 5n, p11)).toBe(5n);
    expect(Fp_add(5n, 0n, p11)).toBe(5n);
  });

  it('should reduce when sum >= p', () => {
    expect(Fp_add(5n, 4n, p7)).toBe(2n); // 9 mod 7 = 2
    expect(Fp_add(6n, 6n, p7)).toBe(5n); // 12 mod 7 = 5
    expect(Fp_add(10n, 10n, p11)).toBe(9n); // 20 mod 11 = 9
  });

  it('should be commutative', () => {
    for (let a = 0n; a < p7; a++) {
      for (let b = 0n; b < p7; b++) {
        expect(Fp_add(a, b, p7)).toBe(Fp_add(b, a, p7));
      }
    }
  });

  it('should be associative', () => {
    expect(Fp_add(Fp_add(2n, 3n, p7), 4n, p7)).toBe(Fp_add(2n, Fp_add(3n, 4n, p7), p7));
  });
});

describe('Fp_sub - Subtraction mod p', () => {
  it('should subtract without underflow', () => {
    expect(Fp_sub(5n, 3n, p7)).toBe(2n);
    expect(Fp_sub(5n, 0n, p11)).toBe(5n);
    expect(Fp_sub(5n, 5n, p11)).toBe(0n);
  });

  it('should reduce when diff < 0', () => {
    expect(Fp_sub(3n, 5n, p7)).toBe(5n); // -2 mod 7 = 5
    expect(Fp_sub(0n, 1n, p7)).toBe(6n); // -1 mod 7 = 6
    expect(Fp_sub(0n, 6n, p7)).toBe(1n); // -6 mod 7 = 1
  });

  it('should satisfy a - b + b = a', () => {
    for (let a = 0n; a < p7; a++) {
      for (let b = 0n; b < p7; b++) {
        expect(Fp_add(Fp_sub(a, b, p7), b, p7)).toBe(a);
      }
    }
  });
});

describe('Fp_neg - Negation mod p', () => {
  it('should negate correctly', () => {
    expect(Fp_neg(0n, p7)).toBe(0n);
    expect(Fp_neg(1n, p7)).toBe(6n);
    expect(Fp_neg(3n, p7)).toBe(4n);
    expect(Fp_neg(6n, p7)).toBe(1n);
  });

  it('should satisfy a + neg(a) = 0', () => {
    for (let a = 0n; a < p11; a++) {
      expect(Fp_add(a, Fp_neg(a, p11), p11)).toBe(0n);
    }
  });

  it('should satisfy neg(neg(a)) = a', () => {
    for (let a = 0n; a < p11; a++) {
      expect(Fp_neg(Fp_neg(a, p11), p11)).toBe(a);
    }
  });
});

/**
 * PARI's Fp_add/Fp_sub/Fp_neg take a fast path when the inputs are already
 * reduced, but always fall back to a full `remii`/`modii` reduction
 * otherwise (pariinl.h:1670-1730), so the result is in [0, p) for arbitrary
 * integer inputs.
 */
describe('Fp_add / Fp_sub / Fp_neg reduce unreduced inputs fully', () => {
  it('matches Fp_red for out-of-range arguments', () => {
    expect(Fp_add(20n, 20n, 11n)).toBe(7n);
    expect(Fp_neg(20n, 11n)).toBe(2n);
    expect(Fp_sub(-5n, 20n, 11n)).toBe(8n);
    expect(Fp_neg(-20n, 11n)).toBe(9n);
    expect(Fp_add(-30n, -30n, 7n)).toBe(3n);
  });

  it('agrees with Fp_red on a dense grid', () => {
    for (const p of [2n, 3n, 5n, 7n, 11n, 13n, 101n]) {
      for (let a = -40n; a <= 40n; a++) {
        expect(Fp_neg(a, p)).toBe(Fp_red(-a, p));
        for (let b = -40n; b <= 40n; b++) {
          expect(Fp_add(a, b, p)).toBe(Fp_red(a + b, p));
          expect(Fp_sub(a, b, p)).toBe(Fp_red(a - b, p));
        }
      }
    }
  });
});

describe('Fp_mul - Multiplication mod p', () => {
  it('should multiply correctly', () => {
    expect(Fp_mul(2n, 3n, p7)).toBe(6n);
    expect(Fp_mul(3n, 3n, p7)).toBe(2n); // 9 mod 7 = 2
    expect(Fp_mul(0n, 5n, p7)).toBe(0n);
    expect(Fp_mul(1n, 5n, p7)).toBe(5n);
  });

  it('should be commutative', () => {
    for (let a = 0n; a < p7; a++) {
      for (let b = 0n; b < p7; b++) {
        expect(Fp_mul(a, b, p7)).toBe(Fp_mul(b, a, p7));
      }
    }
  });

  it('should be associative', () => {
    expect(Fp_mul(Fp_mul(2n, 3n, p7), 4n, p7)).toBe(Fp_mul(2n, Fp_mul(3n, 4n, p7), p7));
  });

  it('should distribute over addition', () => {
    expect(Fp_mul(2n, Fp_add(3n, 4n, p7), p7)).toBe(
      Fp_add(Fp_mul(2n, 3n, p7), Fp_mul(2n, 4n, p7), p7)
    );
  });
});

describe('Fp_sqr - Square mod p', () => {
  it('should square correctly', () => {
    expect(Fp_sqr(0n, p7)).toBe(0n);
    expect(Fp_sqr(1n, p7)).toBe(1n);
    expect(Fp_sqr(2n, p7)).toBe(4n);
    expect(Fp_sqr(3n, p7)).toBe(2n); // 9 mod 7 = 2
    expect(Fp_sqr(4n, p7)).toBe(2n); // 16 mod 7 = 2
  });

  it('should equal a * a', () => {
    for (let a = 0n; a < p11; a++) {
      expect(Fp_sqr(a, p11)).toBe(Fp_mul(a, a, p11));
    }
  });
});

describe('Fp_inv - Modular inverse', () => {
  it('should compute inverses correctly', () => {
    expect(Fp_inv(1n, p7)).toBe(1n);
    expect(Fp_inv(2n, p7)).toBe(4n); // 2 * 4 = 8 = 1 mod 7
    expect(Fp_inv(3n, p7)).toBe(5n); // 3 * 5 = 15 = 1 mod 7
    expect(Fp_inv(4n, p7)).toBe(2n);
    expect(Fp_inv(5n, p7)).toBe(3n);
    expect(Fp_inv(6n, p7)).toBe(6n); // 6 * 6 = 36 = 1 mod 7
  });

  it('should satisfy a * inv(a) = 1', () => {
    for (let a = 1n; a < p11; a++) {
      expect(Fp_mul(a, Fp_inv(a, p11), p11)).toBe(1n);
    }
  });

  it('should satisfy inv(inv(a)) = a', () => {
    for (let a = 1n; a < p11; a++) {
      expect(Fp_inv(Fp_inv(a, p11), p11)).toBe(a);
    }
  });

  it('should throw for non-invertible elements', () => {
    expect(() => Fp_inv(0n, p7)).toThrow();
  });
});

describe('Fp_div - Division mod p', () => {
  it('should divide correctly', () => {
    expect(Fp_div(6n, 2n, p7)).toBe(3n);
    expect(Fp_div(1n, 2n, p7)).toBe(4n); // 1/2 = 4 since 2*4 = 8 = 1 mod 7
    expect(Fp_div(0n, 3n, p7)).toBe(0n);
  });

  it('should satisfy (a/b) * b = a', () => {
    for (let a = 0n; a < p7; a++) {
      for (let b = 1n; b < p7; b++) {
        expect(Fp_mul(Fp_div(a, b, p7), b, p7)).toBe(a);
      }
    }
  });
});

describe('Fp_pow - Modular exponentiation', () => {
  it('should compute powers correctly', () => {
    expect(Fp_pow(2n, 0n, p7)).toBe(1n);
    expect(Fp_pow(2n, 1n, p7)).toBe(2n);
    expect(Fp_pow(2n, 2n, p7)).toBe(4n);
    expect(Fp_pow(2n, 3n, p7)).toBe(1n); // 8 mod 7 = 1
    expect(Fp_pow(3n, 6n, p7)).toBe(1n); // Fermat's little theorem: 3^6 = 1 mod 7
  });

  it('should handle negative exponents', () => {
    expect(Fp_pow(2n, -1n, p7)).toBe(Fp_inv(2n, p7));
    expect(Fp_pow(3n, -2n, p7)).toBe(Fp_inv(Fp_sqr(3n, p7), p7));
  });

  it('should satisfy Fermat little theorem: a^(p-1) = 1 for a != 0', () => {
    for (let a = 1n; a < p11; a++) {
      expect(Fp_pow(a, p11 - 1n, p11)).toBe(1n);
    }
  });

  it('should satisfy a^p = a (mod p) for all a', () => {
    for (let a = 0n; a < p11; a++) {
      expect(Fp_pow(a, p11, p11)).toBe(a);
    }
  });
});

describe('kronecker - Kronecker/Legendre symbol', () => {
  it('should compute Legendre symbol for small values', () => {
    // Quadratic residues mod 7: 1, 2, 4 (since 1^2=1, 2^2=4, 3^2=2)
    expect(kronecker(1n, p7)).toBe(1);
    expect(kronecker(2n, p7)).toBe(1);
    expect(kronecker(3n, p7)).toBe(-1);
    expect(kronecker(4n, p7)).toBe(1);
    expect(kronecker(5n, p7)).toBe(-1);
    expect(kronecker(6n, p7)).toBe(-1);
    expect(kronecker(0n, p7)).toBe(0);
  });

  it('should handle negative values', () => {
    // -1 is a square mod p iff p = 1 mod 4
    expect(kronecker(-1n, p7)).toBe(-1); // 7 = 3 mod 4
    expect(kronecker(-1n, p29)).toBe(1); // 29 = 1 mod 4
    expect(kronecker(-1n, p97)).toBe(1); // 97 = 1 mod 4
  });

  it('should satisfy multiplicativity: (ab|p) = (a|p)(b|p)', () => {
    for (let a = 1n; a < p7; a++) {
      for (let b = 1n; b < p7; b++) {
        expect(kronecker(a * b, p7)).toBe(kronecker(a, p7) * kronecker(b, p7));
      }
    }
  });

  it('should match Euler criterion: (a|p) = a^((p-1)/2) mod p', () => {
    for (let a = 1n; a < p11; a++) {
      const euler = Fp_pow(a, (p11 - 1n) / 2n, p11);
      const expected = euler === 1n ? 1 : -1; // euler is 1 or p-1 (=-1)
      expect(kronecker(a, p11)).toBe(expected);
    }
  });
});

describe('Fp_issquare - Quadratic residue test', () => {
  it('should identify squares correctly', () => {
    // Squares mod 7: 0, 1, 2, 4
    expect(Fp_issquare(0n, p7)).toBe(true);
    expect(Fp_issquare(1n, p7)).toBe(true);
    expect(Fp_issquare(2n, p7)).toBe(true);
    expect(Fp_issquare(3n, p7)).toBe(false);
    expect(Fp_issquare(4n, p7)).toBe(true);
    expect(Fp_issquare(5n, p7)).toBe(false);
    expect(Fp_issquare(6n, p7)).toBe(false);
  });

  it('should be consistent with kronecker', () => {
    for (let a = 1n; a < p11; a++) {
      expect(Fp_issquare(a, p11)).toBe(kronecker(a, p11) === 1);
    }
  });
});

describe('Fp_sqrt - Square root mod p', () => {
  it('should return null for non-squares', () => {
    expect(Fp_sqrt(3n, p7)).toBe(null);
    expect(Fp_sqrt(5n, p7)).toBe(null);
    expect(Fp_sqrt(6n, p7)).toBe(null);
  });

  it('should compute square roots correctly', () => {
    // For squares, verify that sqrt(a)^2 = a
    const squares = [0n, 1n, 2n, 4n];
    for (const a of squares) {
      const r = Fp_sqrt(a, p7);
      expect(r).not.toBe(null);
      if (r !== null) {
        expect(Fp_sqr(r, p7)).toBe(a);
      }
    }
  });

  it('should work for p = 3 mod 4', () => {
    // p = 23 is 3 mod 4, so sqrt(a) = a^((p+1)/4) for quadratic residues
    for (let a = 0n; a < p23; a++) {
      if (Fp_issquare(a, p23)) {
        const r = Fp_sqrt(a, p23);
        expect(r).not.toBe(null);
        if (r !== null) {
          expect(Fp_sqr(r, p23)).toBe(a);
        }
      }
    }
  });

  it('should work for p = 1 mod 4', () => {
    // p = 29 is 1 mod 4, needs Tonelli-Shanks
    for (let a = 0n; a < p29; a++) {
      if (Fp_issquare(a, p29)) {
        const r = Fp_sqrt(a, p29);
        expect(r).not.toBe(null);
        if (r !== null) {
          expect(Fp_sqr(r, p29)).toBe(a);
        }
      }
    }
  });

  it('should work for p = 1 mod 8', () => {
    // p = 97 is 1 mod 8, needs full Tonelli-Shanks
    for (let a = 0n; a < p97; a++) {
      if (Fp_issquare(a, p97)) {
        const r = Fp_sqrt(a, p97);
        expect(r).not.toBe(null);
        if (r !== null) {
          expect(Fp_sqr(r, p97)).toBe(a);
        }
      }
    }
  });

  it('should return the smaller of the two roots', () => {
    // sqrt should return r such that r <= p - r
    for (let a = 1n; a < p23; a++) {
      const r = Fp_sqrt(a, p23);
      if (r !== null && r !== 0n) {
        const other = p23 - r;
        expect(r <= other).toBe(true);
      }
    }
  });

  it('should work for larger primes', () => {
    // Test with p = 1009
    const testValues = [2n, 3n, 4n, 5n, 100n, 500n];
    for (const a of testValues) {
      const r = Fp_sqrt(a, p1009);
      if (r !== null) {
        expect(Fp_sqr(r, p1009)).toBe(Fp_red(a, p1009));
      }
    }
  });
});

describe('Fp_center - Centered reduction', () => {
  it('should center values around 0', () => {
    expect(Fp_center(0n, p7)).toBe(0n);
    expect(Fp_center(1n, p7)).toBe(1n);
    expect(Fp_center(2n, p7)).toBe(2n);
    expect(Fp_center(3n, p7)).toBe(3n);
    expect(Fp_center(4n, p7)).toBe(-3n);
    expect(Fp_center(5n, p7)).toBe(-2n);
    expect(Fp_center(6n, p7)).toBe(-1n);
  });
});

describe('Fp_halve and Fp_double', () => {
  it('should satisfy double(halve(a)) = a for even a', () => {
    for (let a = 0n; a < p11; a++) {
      expect(Fp_double(Fp_halve(a, p11), p11)).toBe(a);
    }
  });

  it('should satisfy halve(double(a)) = a', () => {
    for (let a = 0n; a < p11; a++) {
      expect(Fp_halve(Fp_double(a, p11), p11)).toBe(a);
    }
  });
});

describe('Fp_addmul - Fused add-multiply', () => {
  it('should compute x + y*z correctly', () => {
    expect(Fp_addmul(1n, 2n, 3n, p7)).toBe(Fp_add(1n, Fp_mul(2n, 3n, p7), p7));
    expect(Fp_addmul(0n, 2n, 3n, p7)).toBe(Fp_mul(2n, 3n, p7));
    expect(Fp_addmul(5n, 0n, 3n, p7)).toBe(5n);
    expect(Fp_addmul(5n, 2n, 0n, p7)).toBe(5n);
  });
});

describe('gcd and xgcd', () => {
  it('should compute gcd correctly', () => {
    expect(gcd(12n, 8n)).toBe(4n);
    expect(gcd(17n, 13n)).toBe(1n);
    expect(gcd(100n, 25n)).toBe(25n);
    expect(gcd(0n, 5n)).toBe(5n);
    expect(gcd(5n, 0n)).toBe(5n);
  });

  it('should handle negative inputs', () => {
    expect(gcd(-12n, 8n)).toBe(4n);
    expect(gcd(12n, -8n)).toBe(4n);
    expect(gcd(-12n, -8n)).toBe(4n);
  });

  it('should compute extended gcd correctly', () => {
    const [g, x, y] = xgcd(12n, 8n);
    expect(g).toBe(4n);
    expect(12n * x + 8n * y).toBe(g);
  });

  it('should produce Bezout coefficients', () => {
    const testCases: [bigint, bigint][] = [
      [17n, 13n],
      [100n, 37n],
      [1009n, 97n],
    ];
    for (const [a, b] of testCases) {
      const [g, x, y] = xgcd(a, b);
      expect(a * x + b * y).toBe(g);
    }
  });
});

describe('Field axioms for Fp', () => {
  // Verify that Fp forms a field

  it('should have additive identity (0)', () => {
    for (let a = 0n; a < p7; a++) {
      expect(Fp_add(a, 0n, p7)).toBe(a);
    }
  });

  it('should have multiplicative identity (1)', () => {
    for (let a = 0n; a < p7; a++) {
      expect(Fp_mul(a, 1n, p7)).toBe(a);
    }
  });

  it('should have additive inverses', () => {
    for (let a = 0n; a < p7; a++) {
      expect(Fp_add(a, Fp_neg(a, p7), p7)).toBe(0n);
    }
  });

  it('should have multiplicative inverses for non-zero elements', () => {
    for (let a = 1n; a < p7; a++) {
      expect(Fp_mul(a, Fp_inv(a, p7), p7)).toBe(1n);
    }
  });

  it('should be closed under operations', () => {
    for (let a = 0n; a < p7; a++) {
      for (let b = 0n; b < p7; b++) {
        const sum = Fp_add(a, b, p7);
        const prod = Fp_mul(a, b, p7);
        expect(sum >= 0n && sum < p7).toBe(true);
        expect(prod >= 0n && prod < p7).toBe(true);
      }
    }
  });
});
