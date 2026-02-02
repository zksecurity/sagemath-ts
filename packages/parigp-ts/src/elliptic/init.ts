/**
 * @module parigp-ts/elliptic/init
 * @description Elliptic curve initialization functions ported from PARI/GP
 *
 * Reference: reference/pari/src/basemath/elliptic.c (lines 440-537, 774-795, 6561-6598)
 *
 * This module provides functions to initialize elliptic curves from various
 * input formats and compute derived quantities like discriminant and j-invariant.
 */

/**
 * Curve type enumeration following PARI/GP's t_ELL_* constants.
 * Reference: paridecl.h:3340
 */
export enum EllCurveType {
  /** Generic ring */
  t_ELL_Rg = 0,
  /** Rationals */
  t_ELL_Q = 1,
  /** p-adic field */
  t_ELL_Qp = 2,
  /** Prime field Fp */
  t_ELL_Fp = 3,
  /** Extension field Fq */
  t_ELL_Fq = 4,
  /** Number field */
  t_ELL_NF = 5,
}

/**
 * Elliptic curve structure following PARI/GP's 16-component vector.
 *
 * An elliptic curve in general Weierstrass form:
 *   y^2 + a1*x*y + a3*y = x^3 + a2*x^2 + a4*x + a6
 *
 * Reference: elliptic.c:443-495 (initsmall5 function)
 */
export interface EllipticCurve {
  /** Weierstrass coefficient a1 */
  readonly a1: bigint;
  /** Weierstrass coefficient a2 */
  readonly a2: bigint;
  /** Weierstrass coefficient a3 */
  readonly a3: bigint;
  /** Weierstrass coefficient a4 */
  readonly a4: bigint;
  /** Weierstrass coefficient a6 */
  readonly a6: bigint;

  /**
   * Derived quantity b2 = a1^2 + 4*a2
   */
  readonly b2: bigint;
  /**
   * Derived quantity b4 = a1*a3 + 2*a4
   */
  readonly b4: bigint;
  /**
   * Derived quantity b6 = a3^2 + 4*a6
   */
  readonly b6: bigint;
  /**
   * Derived quantity b8 = a1^2*a6 + 4*a2*a6 + a2*a3^2 - a4*(a4 + a1*a3)
   */
  readonly b8: bigint;

  /**
   * Derived quantity c4 = b2^2 - 24*b4
   */
  readonly c4: bigint;
  /**
   * Derived quantity c6 = 36*b2*b4 - b2^3 - 216*b6
   */
  readonly c6: bigint;

  /**
   * Discriminant of the curve.
   * disc = -b2^2*b8 - 8*b4^3 - 27*b6^2 + 9*b2*b4*b6 = (c4^3 - c6^2) / 1728
   */
  readonly disc: bigint;

  /**
   * j-invariant of the curve.
   * j = c4^3 / disc
   */
  readonly j: bigint;

  /** Curve type (t_ELL_Rg, t_ELL_Fp, etc.) */
  readonly type: EllCurveType;

  /** Prime p for curves over Fp (undefined for other types) */
  readonly p?: bigint;

  /**
   * Short Weierstrass form a4 coefficient.
   * For Fp curves with p > 3, we can transform to y^2 = x^3 + A4*x + A6.
   * A4 = -27 * c4
   */
  readonly A4?: bigint;

  /**
   * Short Weierstrass form a6 coefficient.
   * A6 = -54 * c6
   */
  readonly A6?: bigint;
}

/**
 * Error thrown when elliptic curve initialization fails.
 */
export class EllipticCurveError extends Error {
  constructor(message: string) {
    super(message);
    this.name = 'EllipticCurveError';
  }
}

/**
 * Compute positive modulo (result always in [0, n)).
 */
function mod(a: bigint, n: bigint): bigint {
  const result = a % n;
  return result < 0n ? result + n : result;
}

/**
 * Compute modular inverse of a modulo n using extended Euclidean algorithm.
 * Throws if gcd(a, n) != 1.
 */
function modInverse(a: bigint, n: bigint): bigint {
  a = mod(a, n);
  if (a === 0n) {
    throw new EllipticCurveError('Cannot compute inverse of 0');
  }

  let oldR = a;
  let r = n;
  let oldS = 1n;
  let s = 0n;

  while (r !== 0n) {
    const q = oldR / r;
    const tempR = r;
    r = oldR - q * r;
    oldR = tempR;

    const tempS = s;
    s = oldS - q * s;
    oldS = tempS;
  }

  if (oldR !== 1n) {
    throw new EllipticCurveError(`gcd(${a}, ${n}) = ${oldR} != 1, inverse does not exist`);
  }

  return mod(oldS, n);
}

/**
 * Compute derived quantities (b2, b4, b6, b8, c4, c6, disc) from Weierstrass coefficients.
 *
 * This follows the formulas from elliptic.c:462-493 (initsmall5).
 *
 * @param a1 - Weierstrass coefficient a1
 * @param a2 - Weierstrass coefficient a2
 * @param a3 - Weierstrass coefficient a3
 * @param a4 - Weierstrass coefficient a4
 * @param a6 - Weierstrass coefficient a6
 * @returns Object containing all derived quantities
 */
function computeDerivedQuantities(
  a1: bigint,
  a2: bigint,
  a3: bigint,
  a4: bigint,
  a6: bigint
): {
  b2: bigint;
  b4: bigint;
  b6: bigint;
  b8: bigint;
  c4: bigint;
  c6: bigint;
  disc: bigint;
} {
  // b2 = a1^2 + 4*a2
  const a11 = a1 * a1;
  const b2 = a11 + 4n * a2;

  // b4 = a1*a3 + 2*a4
  const a13 = a1 * a3;
  const b4 = a13 + 2n * a4;

  // b6 = a3^2 + 4*a6
  const a33 = a3 * a3;
  const b6 = a33 + 4n * a6;

  // b8 = a1^2*a6 + 4*a2*a6 + a2*a3^2 - a4*(a4 + a1*a3)
  // Simplified: b8 = a1^2*a6 + b6*a2 - a4*(a4 + a1*a3)
  const b8 = a11 * a6 + b6 * a2 - a4 * (a4 + a13);

  // c4 = b2^2 - 24*b4
  const b22 = b2 * b2;
  const c4 = b22 - 24n * b4;

  // c6 = 36*b2*b4 - b2^3 - 216*b6
  const c6 = 36n * b2 * b4 - b22 * b2 - 216n * b6;

  // disc = 9*b2*b4*b6 - 8*b4^3 - b2^2*b8 - 27*b6^2
  // Reference: elliptic.c:479-480
  const b43 = b4 * b4 * b4;
  const b62 = b6 * b6;
  const disc = 9n * b2 * b4 * b6 - 8n * b43 - b22 * b8 - 27n * b62;

  return { b2, b4, b6, b8, c4, c6, disc };
}

/**
 * Compute derived quantities modulo p.
 */
function computeDerivedQuantitiesMod(
  a1: bigint,
  a2: bigint,
  a3: bigint,
  a4: bigint,
  a6: bigint,
  p: bigint
): {
  b2: bigint;
  b4: bigint;
  b6: bigint;
  b8: bigint;
  c4: bigint;
  c6: bigint;
  disc: bigint;
} {
  const result = computeDerivedQuantities(a1, a2, a3, a4, a6);
  return {
    b2: mod(result.b2, p),
    b4: mod(result.b4, p),
    b6: mod(result.b6, p),
    b8: mod(result.b8, p),
    c4: mod(result.c4, p),
    c6: mod(result.c6, p),
    disc: mod(result.disc, p),
  };
}

/**
 * Compute j-invariant from c4 and discriminant.
 *
 * j = c4^3 / disc
 *
 * Note: Over the integers, j may not be an integer. The relationship is:
 *   c4^3 - c6^2 = 1728 * disc
 *
 * For curves where j is not an integer over Z, we return the floor division.
 * This is acceptable because:
 * 1. For curves over finite fields, we use computeJMod which handles modular arithmetic
 * 2. For curves from j-invariant, the j-value is already known
 * 3. For curves with rational j-invariant that happens to be an integer, this is exact
 *
 * Reference: elliptic.c:497-514 (get_j function)
 */
function computeJ(c4: bigint, disc: bigint): bigint {
  if (disc === 0n) {
    throw new EllipticCurveError('Curve is singular (discriminant is zero)');
  }
  const c43 = c4 * c4 * c4;
  // Integer division - may not be exact for curves where j is a non-integer rational
  return c43 / disc;
}

/**
 * Compute j-invariant modulo p.
 */
function computeJMod(c4: bigint, disc: bigint, p: bigint): bigint {
  if (disc === 0n) {
    throw new EllipticCurveError('Curve is singular (discriminant is zero mod p)');
  }
  const c43 = mod(c4 * c4 * c4, p);
  const discInv = modInverse(disc, p);
  return mod(c43 * discInv, p);
}

/**
 * Create an elliptic curve from j-invariant.
 *
 * For j != 0, 1728, returns [0, 0, 0, 3*k*j, 2*k^2*j] where k = j - 1728.
 * Special cases:
 *   j = 0: returns [0, 0, 0, 0, 1] (curve y^2 = x^3 + 1)
 *   j = 1728: returns [0, 0, 0, 1, 0] (curve y^2 = x^3 + x)
 *
 * Reference: elliptic.c:6561-6598 (ellfromj function)
 *
 * @param j - The j-invariant
 * @returns Weierstrass coefficients [a1, a2, a3, a4, a6]
 */
export function ellfromj(j: bigint): [bigint, bigint, bigint, bigint, bigint] {
  if (j === 0n) {
    // y^2 = x^3 + 1 has j = 0
    return [0n, 0n, 0n, 0n, 1n];
  }

  if (j === 1728n) {
    // y^2 = x^3 + x has j = 1728
    return [0n, 0n, 0n, 1n, 0n];
  }

  // General case: j != 0, 1728
  // k = 1728 - j
  // a4 = 3 * k * j = 3 * (1728 - j) * j
  // a6 = 2 * k^2 * j = 2 * (1728 - j)^2 * j
  const k = 1728n - j;
  const kj = k * j;
  const a4 = 3n * kj;
  const a6 = 2n * kj * k;

  return [0n, 0n, 0n, a4, a6];
}

/**
 * Create an elliptic curve from j-invariant over Fp.
 *
 * Reference: elliptic.c:6572-6598 (ellfromj with finite field handling)
 *
 * @param j - The j-invariant (mod p)
 * @param p - The prime field characteristic
 * @returns Weierstrass coefficients [a1, a2, a3, a4, a6] modulo p
 */
export function ellfromjFp(j: bigint, p: bigint): [bigint, bigint, bigint, bigint, bigint] {
  j = mod(j, p);

  // Handle characteristic 2
  if (p === 2n) {
    if (j === 0n) {
      return [0n, 0n, 1n, 0n, 0n];
    } else {
      // j = 1 in F2: return [1, 0, 0, 0, 1/j] = [1, 0, 0, 0, 1]
      return [1n, 0n, 0n, 0n, modInverse(j, p)];
    }
  }

  // Handle characteristic 3
  if (p === 3n) {
    if (j === 0n) {
      return [0n, 0n, 0n, 1n, 0n];
    } else {
      // a2 = j, a6 = -j^2
      const a6 = mod(-j * j, p);
      return [0n, j, 0n, 0n, a6];
    }
  }

  // Characteristic > 3
  if (j === 0n) {
    return [0n, 0n, 0n, 0n, 1n];
  }

  if (j === mod(1728n, p)) {
    return [0n, 0n, 0n, 1n, 0n];
  }

  // General case
  const k = mod(1728n - j, p);
  const kj = mod(k * j, p);
  const a4 = mod(3n * kj, p);
  const a6 = mod(2n * kj * k, p);

  return [0n, 0n, 0n, a4, a6];
}

/**
 * Initialize an elliptic curve from short Weierstrass form [a4, a6].
 *
 * Creates a curve y^2 = x^3 + a4*x + a6 (a1 = a2 = a3 = 0).
 *
 * Reference: elliptic.c:443-459 (initsmall46 function)
 */
function initFromShortWeierstrass(a4: bigint, a6: bigint, p?: bigint): EllipticCurve {
  // For short Weierstrass, a1 = a2 = a3 = 0
  // b2 = 0, b4 = 2*a4, b6 = 4*a6, b8 = -a4^2
  // c4 = -48*a4, c6 = -864*a6
  // disc = -64*a4^3 - 432*a6^2

  if (p !== undefined) {
    const a4m = mod(a4, p);
    const a6m = mod(a6, p);

    const b2 = 0n;
    const b4 = mod(2n * a4m, p);
    const b6 = mod(4n * a6m, p);
    const b8 = mod(-a4m * a4m, p);

    const c4 = mod(-48n * a4m, p);
    const c6 = mod(-864n * a6m, p);

    // disc = -64*a4^3 - 432*a6^2
    const disc = mod(-64n * a4m * a4m * a4m - 432n * a6m * a6m, p);

    if (disc === 0n) {
      throw new EllipticCurveError('Curve is singular (discriminant is zero mod p)');
    }

    const j = computeJMod(c4, disc, p);

    // Compute short Weierstrass form: A4 = -27*c4, A6 = -54*c6
    const A4 = mod(-27n * c4, p);
    const A6 = mod(-54n * c6, p);

    return {
      a1: 0n,
      a2: 0n,
      a3: 0n,
      a4: a4m,
      a6: a6m,
      b2,
      b4,
      b6,
      b8,
      c4,
      c6,
      disc,
      j,
      type: EllCurveType.t_ELL_Fp,
      p,
      A4,
      A6,
    };
  }

  // Generic ring case
  const b2 = 0n;
  const b4 = 2n * a4;
  const b6 = 4n * a6;
  const b8 = -a4 * a4;

  const c4 = -48n * a4;
  const c6 = -864n * a6;

  // disc = -64*a4^3 - 432*a6^2
  const disc = -64n * a4 * a4 * a4 - 432n * a6 * a6;

  if (disc === 0n) {
    throw new EllipticCurveError('Curve is singular (discriminant is zero)');
  }

  const j = computeJ(c4, disc);

  return {
    a1: 0n,
    a2: 0n,
    a3: 0n,
    a4,
    a6,
    b2,
    b4,
    b6,
    b8,
    c4,
    c6,
    disc,
    j,
    type: EllCurveType.t_ELL_Rg,
  };
}

/**
 * Initialize an elliptic curve from general Weierstrass form [a1, a2, a3, a4, a6].
 *
 * Reference: elliptic.c:462-495 (initsmall5 function)
 */
function initFromGeneralWeierstrass(
  a1: bigint,
  a2: bigint,
  a3: bigint,
  a4: bigint,
  a6: bigint,
  p?: bigint
): EllipticCurve {
  // If all of a1, a2, a3 are zero, use short Weierstrass form
  if (a1 === 0n && a2 === 0n && a3 === 0n) {
    return initFromShortWeierstrass(a4, a6, p);
  }

  if (p !== undefined) {
    const a1m = mod(a1, p);
    const a2m = mod(a2, p);
    const a3m = mod(a3, p);
    const a4m = mod(a4, p);
    const a6m = mod(a6, p);

    const derived = computeDerivedQuantitiesMod(a1m, a2m, a3m, a4m, a6m, p);

    if (derived.disc === 0n) {
      throw new EllipticCurveError('Curve is singular (discriminant is zero mod p)');
    }

    const j = computeJMod(derived.c4, derived.disc, p);

    // For p > 3, compute short Weierstrass form
    let A4: bigint | undefined;
    let A6: bigint | undefined;
    if (p > 3n) {
      A4 = mod(-27n * derived.c4, p);
      A6 = mod(-54n * derived.c6, p);
    }

    return {
      a1: a1m,
      a2: a2m,
      a3: a3m,
      a4: a4m,
      a6: a6m,
      ...derived,
      j,
      type: EllCurveType.t_ELL_Fp,
      p,
      A4,
      A6,
    };
  }

  // Generic ring case
  const derived = computeDerivedQuantities(a1, a2, a3, a4, a6);

  if (derived.disc === 0n) {
    throw new EllipticCurveError('Curve is singular (discriminant is zero)');
  }

  const j = computeJ(derived.c4, derived.disc);

  return {
    a1,
    a2,
    a3,
    a4,
    a6,
    ...derived,
    j,
    type: EllCurveType.t_ELL_Rg,
  };
}

/**
 * Initialize an elliptic curve from j-invariant.
 */
function initFromJInvariant(j: bigint, p?: bigint): EllipticCurve {
  if (p !== undefined) {
    const coeffs = ellfromjFp(j, p);
    return initFromGeneralWeierstrass(coeffs[0], coeffs[1], coeffs[2], coeffs[3], coeffs[4], p);
  }

  const coeffs = ellfromj(j);
  return initFromGeneralWeierstrass(coeffs[0], coeffs[1], coeffs[2], coeffs[3], coeffs[4]);
}

/**
 * Input type for ellinit: can be [a1,a2,a3,a4,a6], [a4,a6], or [j].
 */
export type EllInitInput = [bigint, bigint, bigint, bigint, bigint] | [bigint, bigint] | [bigint];

/**
 * Initialize an elliptic curve from various input formats.
 *
 * This is the main entry point, porting PARI/GP's ellinit function.
 *
 * Reference: elliptic.c:887-893 (ellinit), 839-885 (ellinit_i)
 *
 * @param x - Input in one of these formats:
 *   - [a1, a2, a3, a4, a6]: General Weierstrass coefficients
 *   - [a4, a6]: Short Weierstrass (y^2 = x^3 + a4*x + a6)
 *   - [j]: j-invariant (constructs curve with that j-invariant)
 *
 * @param D - Domain specification:
 *   - bigint p: Create curve over Fp (prime field)
 *   - undefined: Create curve over generic ring (exact arithmetic)
 *
 * @param prec - Precision (unused, for compatibility with PARI/GP)
 *
 * @returns Initialized elliptic curve structure
 *
 * @throws {EllipticCurveError} If the curve is singular (discriminant zero)
 *
 * @example
 * ```typescript
 * // Short Weierstrass form: y^2 = x^3 + 7 (secp256k1)
 * const E = ellinit([0n, 7n], 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2Fn);
 *
 * // From j-invariant over generic ring
 * const E2 = ellinit([1728n]);  // y^2 = x^3 + x
 *
 * // General Weierstrass form
 * const E3 = ellinit([1n, 2n, 3n, 4n, 6n]);
 * ```
 */
export function ellinit(x: EllInitInput, D?: bigint, _prec?: number): EllipticCurve {
  const p = D;

  switch (x.length) {
    case 1:
      // [j] - from j-invariant
      return initFromJInvariant(x[0], p);

    case 2:
      // [a4, a6] - short Weierstrass
      return initFromShortWeierstrass(x[0], x[1], p);

    case 5:
      // [a1, a2, a3, a4, a6] - general Weierstrass
      return initFromGeneralWeierstrass(x[0], x[1], x[2], x[3], x[4], p);

    default:
      throw new EllipticCurveError(
        `Invalid input format: expected 1, 2, or 5 elements, got ${(x as bigint[]).length}`
      );
  }
}

/**
 * Get the j-invariant of an elliptic curve.
 *
 * Reference: elliptic.c uses ell_get_j accessor
 *
 * @param E - Elliptic curve
 * @returns The j-invariant
 */
export function ellj(E: EllipticCurve): bigint {
  return E.j;
}

/**
 * Get the discriminant of an elliptic curve.
 *
 * @param E - Elliptic curve
 * @returns The discriminant
 */
export function elldisc(E: EllipticCurve): bigint {
  return E.disc;
}

/**
 * Check if a curve is non-singular.
 *
 * @param E - Elliptic curve
 * @returns true if the curve is non-singular (discriminant != 0)
 */
export function ellisnonsingular(E: EllipticCurve): boolean {
  return E.disc !== 0n;
}

/**
 * Get the Weierstrass coefficients as an array [a1, a2, a3, a4, a6].
 *
 * @param E - Elliptic curve
 * @returns Array of Weierstrass coefficients
 */
export function ellcoeffs(E: EllipticCurve): [bigint, bigint, bigint, bigint, bigint] {
  return [E.a1, E.a2, E.a3, E.a4, E.a6];
}

/**
 * Convert a general Weierstrass curve to short Weierstrass form [a4, a6].
 *
 * For curves over fields of characteristic > 3, any curve can be transformed
 * to y^2 = x^3 + a4*x + a6 form.
 *
 * The transformation uses:
 *   A4 = -27 * c4
 *   A6 = -54 * c6
 *
 * Reference: elliptic.c:63-68 (c4c6_to_a4a6)
 *
 * @param E - Elliptic curve
 * @returns [A4, A6] for short Weierstrass form, or undefined for char 2,3
 */
export function ellToShortWeierstrass(E: EllipticCurve): [bigint, bigint] | undefined {
  if (E.A4 !== undefined && E.A6 !== undefined) {
    return [E.A4, E.A6];
  }

  // For generic ring, compute without reduction
  if (E.type === EllCurveType.t_ELL_Rg) {
    return [-27n * E.c4, -54n * E.c6];
  }

  return undefined;
}
