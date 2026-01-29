/**
 * @module sage/rings/finite_rings/conway_polynomials
 * @description Database of Conway polynomials for standardized irreducible polynomials
 *
 * Conway polynomials provide a canonical choice of irreducible polynomial
 * for constructing finite field extensions. They have special compatibility
 * properties that make them useful for embedding smaller fields into larger ones.
 *
 * Reference: https://www.math.rwth-aachen.de/~Frank.Luebeck/data/ConwayPol/index.html
 *
 * Format: CONWAY_POLYNOMIALS[p][n] = [c_0, c_1, ..., c_{n-1}]
 * represents the polynomial x^n + c_{n-1}*x^{n-1} + ... + c_1*x + c_0
 * (coefficients are stored in increasing degree order, monic polynomial so leading coeff is 1)
 */

/**
 * Conway polynomials indexed by [characteristic][degree].
 * Coefficients are stored in increasing degree order (constant term first).
 * The polynomial is monic, so the leading coefficient (x^n) is implicitly 1.
 */
export const CONWAY_POLYNOMIALS: Record<number, Record<number, number[]>> = {
  // GF(2^n) - Binary field extensions
  2: {
    // x^2 + x + 1
    2: [1, 1],
    // x^3 + x + 1
    3: [1, 1, 0],
    // x^4 + x + 1
    4: [1, 1, 0, 0],
    // x^5 + x^2 + 1
    5: [1, 0, 1, 0, 0],
    // x^6 + x^4 + x^3 + x + 1
    6: [1, 1, 0, 1, 1, 0],
    // x^7 + x + 1
    7: [1, 1, 0, 0, 0, 0, 0],
    // x^8 + x^4 + x^3 + x^2 + 1
    8: [1, 0, 1, 1, 1, 0, 0, 0],
    // x^9 + x^4 + 1
    9: [1, 0, 0, 0, 1, 0, 0, 0, 0],
    // x^10 + x^6 + x^5 + x^3 + x^2 + x + 1
    10: [1, 1, 1, 1, 0, 1, 1, 0, 0, 0],
    // x^11 + x^2 + 1
    11: [1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
    // x^12 + x^7 + x^6 + x^5 + x^3 + x + 1
    12: [1, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0],
    // x^13 + x^4 + x^3 + x + 1
    13: [1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0],
    // x^14 + x^12 + x^11 + x + 1
    14: [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0],
    // x^15 + x^5 + x^4 + x^2 + 1
    15: [1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    // x^16 + x^5 + x^3 + x^2 + 1
    16: [1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    // x^17 + x^3 + 1
    17: [1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    // x^18 + x^12 + x^10 + x + 1
    18: [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0],
    // x^19 + x^5 + x^2 + x + 1
    19: [1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    // x^20 + x^10 + x^9 + x^7 + x^6 + x^5 + x^4 + x + 1
    20: [1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    // x^32 + x^7 + x^5 + x^3 + x^2 + x + 1
    32: [
      1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0,
    ],
    // x^64 + x^4 + x^3 + x + 1
    64: new Array(64)
      .fill(0)
      .map((_, i) => (i === 0 ? 1 : i === 1 ? 1 : i === 3 ? 1 : i === 4 ? 1 : 0)),
    // x^128 + x^7 + x^2 + x + 1
    128: new Array(128)
      .fill(0)
      .map((_, i) => (i === 0 ? 1 : i === 1 ? 1 : i === 2 ? 1 : i === 7 ? 1 : 0)),
  },

  // GF(3^n) - Ternary field extensions
  3: {
    // x^2 + 2x + 2
    2: [2, 2],
    // x^3 + 2x + 1
    3: [1, 2, 0],
    // x^4 + 2x^3 + 2
    4: [2, 0, 0, 2],
    // x^5 + 2x + 1
    5: [1, 2, 0, 0, 0],
    // x^6 + 2x^4 + x^2 + 2x + 2
    6: [2, 2, 1, 0, 2, 0],
    // x^7 + 2x^2 + 1
    7: [1, 0, 2, 0, 0, 0, 0],
    // x^8 + 2x^5 + x^4 + 2x^2 + 2
    8: [2, 0, 2, 0, 1, 2, 0, 0],
  },

  // GF(5^n)
  5: {
    // x^2 + 4x + 2
    2: [2, 4],
    // x^3 + 3x + 3
    3: [3, 3, 0],
    // x^4 + 4x^2 + 4x + 2
    4: [2, 4, 4, 0],
    // x^5 + 4x + 2
    5: [2, 4, 0, 0, 0],
    // x^6 + 4x^4 + 3x^3 + 2x^2 + 3
    6: [3, 0, 2, 3, 4, 0],
  },

  // GF(7^n)
  7: {
    // x^2 + 6x + 3
    2: [3, 6],
    // x^3 + 6x^2 + 4
    3: [4, 0, 6],
    // x^4 + 5x^2 + 4x + 3
    4: [3, 4, 5, 0],
    // x^5 + 4x^2 + 6x + 3
    5: [3, 6, 4, 0, 0],
  },

  // GF(11^n)
  11: {
    // x^2 + 7x + 2
    2: [2, 7],
    // x^3 + 2x + 9
    3: [9, 2, 0],
    // x^4 + 8x^2 + 10x + 2
    4: [2, 10, 8, 0],
  },

  // GF(13^n)
  13: {
    // x^2 + 12x + 2
    2: [2, 12],
    // x^3 + 2x + 11
    3: [11, 2, 0],
    // x^4 + 12x^2 + 6x + 2
    4: [2, 6, 12, 0],
  },

  // GF(17^n)
  17: {
    // x^2 + 14x + 3
    2: [3, 14],
    // x^3 + 16x^2 + 3
    3: [3, 0, 16],
  },

  // GF(19^n)
  19: {
    // x^2 + 18x + 2
    2: [2, 18],
    // x^3 + 4x + 17
    3: [17, 4, 0],
  },

  // GF(23^n)
  23: {
    // x^2 + 21x + 5
    2: [5, 21],
    // x^3 + 2x + 21
    3: [21, 2, 0],
  },

  // GF(29^n)
  29: {
    // x^2 + 27x + 2
    2: [2, 27],
    // x^3 + 2x + 27
    3: [27, 2, 0],
  },

  // GF(31^n)
  31: {
    // x^2 + 28x + 3
    2: [3, 28],
    // x^3 + 28x^2 + 3
    3: [3, 0, 28],
  },
};

/**
 * Get the Conway polynomial for GF(p^n).
 *
 * @param p - Prime characteristic
 * @param n - Extension degree
 * @returns Array of coefficients [c_0, c_1, ..., c_{n-1}] for x^n + c_{n-1}*x^{n-1} + ... + c_0
 * @throws Error if Conway polynomial is not in the database
 */
export function conway_polynomial(p: number, n: number): number[] {
  const pPolys = CONWAY_POLYNOMIALS[p];
  if (!pPolys) {
    throw new Error(`No Conway polynomials in database for characteristic ${p}`);
  }

  const coeffs = pPolys[n];
  if (!coeffs) {
    throw new Error(`No Conway polynomial in database for GF(${p}^${n})`);
  }

  return coeffs;
}

/**
 * Check if a Conway polynomial exists in the database.
 *
 * @param p - Prime characteristic
 * @param n - Extension degree
 * @returns true if the Conway polynomial is available
 */
export function has_conway_polynomial(p: number, n: number): boolean {
  const pPolys = CONWAY_POLYNOMIALS[p];
  if (!pPolys) {
    return false;
  }
  return n in pPolys;
}

/**
 * List available degrees for a given characteristic.
 *
 * @param p - Prime characteristic
 * @returns Array of available extension degrees
 */
export function available_degrees(p: number): number[] {
  const pPolys = CONWAY_POLYNOMIALS[p];
  if (!pPolys) {
    return [];
  }
  return Object.keys(pPolys)
    .map(Number)
    .sort((a, b) => a - b);
}

/**
 * List available characteristics in the database.
 *
 * @returns Array of available prime characteristics
 */
export function available_characteristics(): number[] {
  return Object.keys(CONWAY_POLYNOMIALS)
    .map(Number)
    .sort((a, b) => a - b);
}
