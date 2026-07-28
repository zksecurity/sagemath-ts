/**
 * @module sage/rings/finite_rings/conway_polynomials
 * @description Database of Conway polynomials for standardized irreducible polynomials
 *
 * Conway polynomials provide a canonical choice of irreducible polynomial
 * for constructing finite field extensions. They have special compatibility
 * properties that make them useful for embedding smaller fields into larger ones.
 *
 * Every entry is monic, irreducible **and** primitive (the class of `x` generates
 * the multiplicative group of GF(p^n)), and is the lexicographically least such
 * polynomial compatible with all proper subfields (Frank Luebeck's normalisation).
 *
 * Data source: the values are decoded from the table FLINT ships in
 * `reference/flint/src/nmod_poly/conway_polynomial_data.c`
 * (`__nmod_poly_cp_primes0` / `_cp_degrees0` / `_numntcoeffs0` / `_ntcoeffs0`,
 * read by `reference/flint/src/nmod_poly/conway.c::conway_polynomial_lt_260`),
 * which is itself generated from Frank Luebeck's `CPimport.txt`, the same source
 * SageMath's `conway_polynomials` package uses.
 *
 * Reference: https://www.math.rwth-aachen.de/~Frank.Luebeck/data/ConwayPol/index.html
 *
 * Format: CONWAY_POLYNOMIALS[p][n] = [c_0, c_1, ..., c_{n-1}]
 * represents the polynomial x^n + c_{n-1}*x^{n-1} + ... + c_1*x + c_0
 * (coefficients are stored in increasing degree order, monic polynomial so leading coeff is 1)
 */

/**
 * Compact storage, mirroring FLINT's encoding: for each (p, n) only the
 * coefficients up to the highest nonzero one below x^n are listed; every
 * omitted coefficient c_i (i < n) is zero and the leading coefficient is 1.
 *
 * See `conway.c::conway_polynomial_lt_260`, which likewise writes
 * `op[deg] = 1`, zeroes `op[1..deg-1]` and then copies the stored
 * `num_nontrivialcoeffs[kx]` low coefficients over the top.
 */
const CONWAY_LOW_COEFFICIENTS: Record<number, Record<number, number[]>> = {
  2: {
    2: [1, 1],
    3: [1, 1],
    4: [1, 1],
    5: [1, 0, 1],
    6: [1, 1, 0, 1, 1],
    7: [1, 1],
    8: [1, 0, 1, 1, 1],
    9: [1, 0, 0, 0, 1],
    10: [1, 1, 1, 1, 0, 1, 1],
    11: [1, 0, 1],
    12: [1, 1, 0, 1, 0, 1, 1, 1],
    13: [1, 1, 0, 1, 1],
    14: [1, 0, 0, 1, 0, 1, 0, 1],
    15: [1, 0, 1, 0, 1, 1],
    16: [1, 0, 1, 1, 0, 1],
    17: [1, 0, 0, 1],
    18: [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1],
    19: [1, 1, 1, 0, 0, 1],
    20: [1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1],
    21: [1, 0, 1, 0, 0, 1, 1],
    22: [1, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1],
    23: [1, 0, 0, 0, 0, 1],
    24: [1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 1, 1],
    25: [1, 0, 1, 0, 0, 0, 1, 0, 1],
    26: [1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1],
    27: [1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1],
    28: [1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1],
    29: [1, 0, 1],
    30: [1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1],
    31: [1, 0, 0, 1],
    32: [1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1],
    33: [1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1],
    34: [1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1],
    35: [1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1],
    36: [1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1],
    37: [1, 1, 1, 1, 1, 1],
    38: [1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 1],
    39: [1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1],
    40: [1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1],
    41: [1, 0, 0, 1],
    42: [1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1],
    43: [1, 0, 0, 1, 1, 0, 1],
    44: [1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1],
    45: [1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1],
    46: [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 1],
    47: [1, 0, 0, 0, 0, 1],
    48: [1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1],
    49: [1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1],
    50: [1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1],
    51: [1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1],
    52: [1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1],
    53: [1, 1, 1, 0, 0, 0, 1],
    54: [1, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1],
    55: [1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1],
    56: [1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1],
    57: [1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1],
    58: [1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1],
    59: [1, 1, 0, 1, 1, 1, 1],
    60: [1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1],
    61: [1, 1, 1, 0, 0, 1],
    62: [1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1],
    63: [1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1],
    64: [1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 1],
  },
  3: {
    2: [2, 2],
    3: [1, 2],
    4: [2, 0, 0, 2],
    5: [1, 2],
    6: [2, 2, 1, 0, 2],
    7: [1, 0, 2],
    8: [2, 2, 2, 0, 1, 2],
    9: [1, 1, 2, 2],
    10: [2, 1, 0, 0, 2, 2, 2],
    11: [1, 0, 2],
    12: [2, 0, 1, 0, 1, 1, 1],
    13: [1, 2],
    14: [2, 0, 1, 2, 0, 1, 2, 1, 1, 2],
    15: [1, 1, 2, 0, 0, 1, 0, 0, 2],
    16: [2, 1, 2, 2, 2, 0, 2, 2],
    17: [1, 2],
    18: [2, 0, 2, 0, 2, 1, 2, 0, 2, 0, 1],
    19: [1, 0, 2],
    20: [2, 1, 0, 2, 2, 2, 0, 0, 1, 1, 1, 1, 0, 2],
    21: [1, 2, 0, 2, 0, 1, 2, 0, 2, 0, 2],
    22: [2, 2, 0, 1, 0, 1, 1, 1, 2, 2, 1, 2],
    23: [1, 1, 0, 1],
    24: [2, 2, 0, 2, 2, 0, 2, 0, 2, 0, 0, 2, 0, 0, 1],
  },
  5: {
    2: [2, 4],
    3: [3, 3],
    4: [2, 4, 4],
    5: [3, 4],
    6: [2, 0, 1, 4, 1],
    7: [3, 3],
    8: [2, 4, 3, 0, 1],
    9: [3, 1, 0, 2],
    10: [2, 1, 4, 2, 3, 3],
    11: [3, 3],
    12: [2, 2, 3, 4, 4, 0, 1, 1],
    13: [3, 3, 4],
    14: [2, 1, 0, 3, 2, 4, 4, 0, 1],
    15: [3, 4, 3, 3, 0, 2],
    16: [2, 1, 4, 4, 2, 4, 4, 4, 1],
    17: [3, 2, 3],
    18: [2, 0, 2, 2, 0, 1, 2, 0, 2, 1, 1, 1, 1],
  },
  7: {
    2: [3, 6],
    3: [4, 0, 6],
    4: [3, 4, 5],
    5: [4, 1],
    6: [3, 6, 4, 5, 1],
    7: [4, 6],
    8: [3, 2, 6, 4],
    9: [4, 6, 0, 1, 6],
    10: [3, 3, 2, 1, 4, 1, 1],
    11: [4, 1],
    12: [3, 0, 5, 0, 4, 2, 3, 5, 2],
    13: [4, 0, 6],
    14: [3, 6, 3, 0, 2, 6, 0, 5],
  },
  11: {
    2: [2, 7],
    3: [9, 2],
    4: [2, 10, 8],
    5: [9, 0, 10],
    6: [2, 7, 6, 4, 3],
    7: [9, 4],
    8: [2, 7, 1, 7, 7],
    9: [9, 8, 9],
    10: [2, 6, 6, 10, 8, 7],
    11: [9, 10],
    12: [2, 5, 6, 5, 5, 2, 4, 1, 1],
  },
  13: {
    2: [2, 12],
    3: [11, 2],
    4: [2, 12, 3],
    5: [11, 4],
    6: [2, 11, 11, 10],
    7: [11, 3],
    8: [2, 3, 2, 12, 8],
    9: [11, 12, 12, 8, 12],
    10: [2, 1, 1, 8, 5, 7],
    11: [11, 3],
    12: [2, 4, 1, 1, 3, 11, 8, 5, 1],
  },
  17: {
    2: [3, 16],
    3: [14, 1],
    4: [3, 10, 7],
    5: [14, 1],
    6: [3, 3, 10, 0, 2],
    7: [14, 12],
    8: [3, 6, 0, 12, 11],
    9: [14, 8, 7],
    10: [3, 12, 9, 5, 6, 13],
  },
  19: {
    2: [2, 18],
    3: [17, 4],
    4: [2, 11, 2],
    5: [17, 5],
    6: [2, 6, 17, 17],
    7: [17, 6],
    8: [2, 3, 10, 12, 1],
    9: [17, 16, 14, 11],
    10: [2, 4, 3, 17, 13, 18],
  },
  23: {
    2: [5, 21],
    3: [18, 2],
    4: [5, 19, 3],
    5: [18, 3],
    6: [5, 1, 9, 9, 1],
    7: [18, 21],
    8: [5, 3, 5, 20, 3],
    9: [18, 9, 8, 3],
    10: [5, 1, 6, 15, 5, 17],
  },
  29: {
    2: [2, 24],
    3: [27, 2],
    4: [2, 15, 2],
    5: [27, 3],
    6: [2, 13, 17, 25, 1],
    7: [27, 2],
    8: [2, 23, 26, 24, 3],
    9: [27, 22, 22, 4],
    10: [2, 22, 2, 17, 8, 25, 1],
  },
  31: {
    2: [3, 29],
    3: [28, 1],
    4: [3, 16, 3],
    5: [28, 7],
    6: [3, 8, 16, 19],
    7: [28, 1],
    8: [3, 24, 12, 25],
    9: [28, 29, 20, 4],
    10: [3, 13, 13, 13, 26, 30],
  },
};

/**
 * Conway polynomials indexed by [characteristic][degree].
 * Coefficients are stored in increasing degree order (constant term first),
 * zero-padded to exactly `n` entries.
 * The polynomial is monic, so the leading coefficient (x^n) is implicitly 1.
 */
export const CONWAY_POLYNOMIALS: Record<number, Record<number, number[]>> = (() => {
  const table: Record<number, Record<number, number[]>> = {};
  for (const pKey of Object.keys(CONWAY_LOW_COEFFICIENTS)) {
    const p = Number(pKey);
    const byDegree: Record<number, number[]> = {};
    const low = CONWAY_LOW_COEFFICIENTS[p]!;
    for (const nKey of Object.keys(low)) {
      const n = Number(nKey);
      const coeffs = new Array<number>(n).fill(0);
      const stored = low[n]!;
      for (let i = 0; i < stored.length; i++) {
        coeffs[i] = stored[i]!;
      }
      byDegree[n] = coeffs;
    }
    table[p] = byDegree;
  }
  return table;
})();

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
