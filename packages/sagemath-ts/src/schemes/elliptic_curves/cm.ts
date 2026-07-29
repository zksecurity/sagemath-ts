/**
 * @module sage/schemes/elliptic_curves/cm
 * @description Complex multiplication for elliptic curves
 *
 * Port of: sage/schemes/elliptic_curves/cm.py
 *
 * This module implements functions related to complex multiplication
 * for elliptic curves, including Hilbert class polynomials and
 * CM j-invariants.
 */

import { INV_J, polclass0, polmodular_db_init } from '@sagemath-ts/parigp-ts';
import { NotImplementedError, ValueError } from '../../errors.js';
import { Integer, ZZ } from '../../rings/integer_ring.js';
import type { Polynomial } from '../../rings/polynomial/polynomial_element.js';
import { PolynomialRing } from '../../rings/polynomial/polynomial_ring.js';

/**
 * Table from Mark Watkins paper "Class numbers of imaginary quadratic fields".
 *
 * The keys are integers 1--100 and the value for each h is [|D|, n]
 * where |D| is the largest absolute discriminant of an imaginary
 * quadratic field with class number h, and n is the number of such
 * fields. These are all unconditional (not dependent on GRH).
 *
 * @see Reference: sage/schemes/elliptic_curves/cm.py:watkins_table
 */
const watkins_table: Map<number, [bigint, number]> = new Map([
  [1, [163n, 9]],
  [2, [427n, 18]],
  [3, [907n, 16]],
  [4, [1555n, 54]],
  [5, [2683n, 25]],
  [6, [3763n, 51]],
  [7, [5923n, 31]],
  [8, [6307n, 131]],
  [9, [10627n, 34]],
  [10, [13843n, 87]],
  [11, [15667n, 41]],
  [12, [17803n, 206]],
  [13, [20563n, 37]],
  [14, [30067n, 95]],
  [15, [34483n, 68]],
  [16, [31243n, 322]],
  [17, [37123n, 45]],
  [18, [48427n, 150]],
  [19, [38707n, 47]],
  [20, [58507n, 350]],
  [21, [61483n, 85]],
  [22, [85507n, 139]],
  [23, [90787n, 68]],
  [24, [111763n, 511]],
  [25, [93307n, 95]],
  [26, [103027n, 190]],
  [27, [103387n, 93]],
  [28, [126043n, 457]],
  [29, [166147n, 83]],
  [30, [134467n, 255]],
  [31, [133387n, 73]],
  [32, [164803n, 708]],
  [33, [222643n, 101]],
  [34, [189883n, 219]],
  [35, [210907n, 103]],
  [36, [217627n, 668]],
  [37, [158923n, 85]],
  [38, [289963n, 237]],
  [39, [253507n, 115]],
  [40, [260947n, 912]],
  [41, [296587n, 109]],
  [42, [280267n, 339]],
  [43, [300787n, 106]],
  [44, [319867n, 691]],
  [45, [308323n, 154]],
  [46, [462883n, 268]],
  [47, [375523n, 107]],
  [48, [335203n, 1365]],
  [49, [393187n, 132]],
  [50, [389467n, 345]],
  [51, [546067n, 159]],
  [52, [439147n, 770]],
  [53, [425107n, 114]],
  [54, [532123n, 427]],
  [55, [452083n, 163]],
  [56, [494323n, 1205]],
  [57, [615883n, 179]],
  [58, [586987n, 291]],
  [59, [474307n, 128]],
  [60, [662803n, 1302]],
  [61, [606643n, 132]],
  [62, [647707n, 323]],
  [63, [991027n, 216]],
  [64, [693067n, 1672]],
  [65, [703123n, 164]],
  [66, [958483n, 530]],
  [67, [652723n, 120]],
  [68, [819163n, 976]],
  [69, [888427n, 209]],
  [70, [811507n, 560]],
  [71, [909547n, 150]],
  [72, [947923n, 1930]],
  [73, [886867n, 119]],
  [74, [951043n, 407]],
  [75, [916507n, 237]],
  [76, [1086187n, 1075]],
  [77, [1242763n, 216]],
  [78, [1004347n, 561]],
  [79, [1333963n, 175]],
  [80, [1165483n, 2277]],
  [81, [1030723n, 228]],
  [82, [1446547n, 402]],
  [83, [1074907n, 150]],
  [84, [1225387n, 1715]],
  [85, [1285747n, 221]],
  [86, [1534723n, 472]],
  [87, [1261747n, 222]],
  [88, [1265587n, 1905]],
  [89, [1429387n, 192]],
  [90, [1548523n, 801]],
  [91, [1391083n, 214]],
  [92, [1452067n, 1248]],
  [93, [1475203n, 262]],
  [94, [1587763n, 509]],
  [95, [1659067n, 241]],
  [96, [1684027n, 3283]],
  [97, [1842523n, 185]],
  [98, [2383747n, 580]],
  [99, [1480627n, 289]],
  [100, [1856563n, 1736]],
]);

/**
 * Table from Janis Klaise [Klaise2012].
 *
 * The keys are integers 1--100 and the value for each h is [|D|, n]
 * where |D| is the largest discriminant of an imaginary quadratic
 * order with class number h, and n is the number of such orders.
 * These are all unconditional (not dependent on GRH).
 *
 * @see Reference: sage/schemes/elliptic_curves/cm.py:klaise_table
 */
const klaise_table: Map<number, [bigint, number]> = new Map([
  [1, [163n, 13]],
  [2, [427n, 29]],
  [3, [907n, 25]],
  [4, [1555n, 84]],
  [5, [2683n, 29]],
  [6, [4075n, 101]],
  [7, [5923n, 38]],
  [8, [7987n, 208]],
  [9, [10627n, 55]],
  [10, [13843n, 123]],
  [11, [15667n, 46]],
  [12, [19723n, 379]],
  [13, [20563n, 43]],
  [14, [30067n, 134]],
  [15, [34483n, 95]],
  [16, [35275n, 531]],
  [17, [37123n, 50]],
  [18, [48427n, 291]],
  [19, [38707n, 59]],
  [20, [58843n, 502]],
  [21, [61483n, 118]],
  [22, [85507n, 184]],
  [23, [90787n, 78]],
  [24, [111763n, 1042]],
  [25, [93307n, 101]],
  [26, [103027n, 227]],
  [27, [103387n, 136]],
  [28, [126043n, 623]],
  [29, [166147n, 94]],
  [30, [137083n, 473]],
  [31, [133387n, 83]],
  [32, [164803n, 1231]],
  [33, [222643n, 158]],
  [34, [189883n, 262]],
  [35, [210907n, 111]],
  [36, [217627n, 1306]],
  [37, [158923n, 96]],
  [38, [289963n, 284]],
  [39, [253507n, 162]],
  [40, [274003n, 1418]],
  [41, [296587n, 125]],
  [42, [301387n, 596]],
  [43, [300787n, 123]],
  [44, [319867n, 911]],
  [45, [308323n, 231]],
  [46, [462883n, 330]],
  [47, [375523n, 117]],
  [48, [335203n, 2895]],
  [49, [393187n, 146]],
  [50, [389467n, 445]],
  [51, [546067n, 217]],
  [52, [457867n, 1006]],
  [53, [425107n, 130]],
  [54, [532123n, 812]],
  [55, [452083n, 177]],
  [56, [494323n, 1812]],
  [57, [615883n, 237]],
  [58, [586987n, 361]],
  [59, [474307n, 144]],
  [60, [662803n, 2361]],
  [61, [606643n, 149]],
  [62, [647707n, 386]],
  [63, [991027n, 311]],
  [64, [693067n, 2919]],
  [65, [703123n, 192]],
  [66, [958483n, 861]],
  [67, [652723n, 145]],
  [68, [819163n, 1228]],
  [69, [888427n, 292]],
  [70, [821683n, 704]],
  [71, [909547n, 176]],
  [72, [947923n, 4059]],
  [73, [886867n, 137]],
  [74, [951043n, 474]],
  [75, [916507n, 353]],
  [76, [1086187n, 1384]],
  [77, [1242763n, 236]],
  [78, [1004347n, 925]],
  [79, [1333963n, 200]],
  [80, [1165483n, 3856]],
  [81, [1030723n, 339]],
  [82, [1446547n, 487]],
  [83, [1074907n, 174]],
  [84, [1225387n, 2998]],
  [85, [1285747n, 246]],
  [86, [1534723n, 555]],
  [87, [1261747n, 313]],
  [88, [1265587n, 2771]],
  [89, [1429387n, 206]],
  [90, [1548523n, 1516]],
  [91, [1391083n, 249]],
  [92, [1452067n, 1591]],
  [93, [1475203n, 354]],
  [94, [1587763n, 600]],
  [95, [1659067n, 273]],
  [96, [1684027n, 7276]],
  [97, [1842523n, 208]],
  [98, [2383747n, 710]],
  [99, [1480627n, 396]],
  [100, [1856563n, 2311]],
]);

/**
 * Dictionary of CM orders by class number.
 * Keys are class numbers, values are lists of [D, f] pairs.
 * Initialized with h=1 data.
 */
const hDf_dict: Map<number, Array<[bigint, bigint]>> = new Map([
  [
    1,
    [
      [-3n, 1n],
      [-3n, 2n],
      [-3n, 3n],
      [-4n, 1n],
      [-4n, 2n],
      [-7n, 1n],
      [-7n, 2n],
      [-8n, 1n],
      [-11n, 1n],
      [-19n, 1n],
      [-43n, 1n],
      [-67n, 1n],
      [-163n, 1n],
    ],
  ],
]);

/**
 * CM j-invariants over Q with their discriminants.
 * Format: [discriminant, conductor, j-invariant]
 */
const cm_j_invariants_QQ: Array<[bigint, bigint, bigint]> = [
  [-3n, 3n, -12288000n],
  [-3n, 2n, 54000n],
  [-3n, 1n, 0n],
  [-4n, 2n, 287496n],
  [-4n, 1n, 1728n],
  [-7n, 2n, 16581375n],
  [-7n, 1n, -3375n],
  [-8n, 1n, 8000n],
  [-11n, 1n, -32768n],
  [-19n, 1n, -884736n],
  [-43n, 1n, -884736000n],
  [-67n, 1n, -147197952000n],
  [-163n, 1n, -262537412640768000n],
];

/**
 * Compute the Kronecker symbol (D/n).
 *
 * @param D - discriminant
 * @param n - positive integer
 * @returns -1, 0, or 1
 */
function kronecker_symbol(D: bigint, n: bigint): bigint {
  if (n === 0n) return D === 1n || D === -1n ? 1n : 0n;
  if (n < 0n) {
    return D < 0n ? -kronecker_symbol(D, -n) : kronecker_symbol(D, -n);
  }

  // Handle n = 1
  if (n === 1n) return 1n;

  // Factor out powers of 2
  let v = 0n;
  let nTemp = n;
  while ((nTemp & 1n) === 0n) {
    v++;
    nTemp >>= 1n;
  }

  // Compute (D/2^v) * (D/odd part)
  let result = 1n;

  // (D/2) depends on D mod 8
  if (v > 0n) {
    const Dmod8 = ((D % 8n) + 8n) % 8n;
    if (Dmod8 === 1n || Dmod8 === 7n) {
      // (D/2) = 1
    } else if (Dmod8 === 3n || Dmod8 === 5n) {
      if ((v & 1n) === 1n) result = -result;
    } else {
      result = 0n;
    }
  }

  if (result === 0n) return 0n;

  // Now compute (D/nTemp) where nTemp is odd
  // Using quadratic reciprocity
  let a = ((D % nTemp) + nTemp) % nTemp;
  let b = nTemp;

  while (a !== 0n) {
    while ((a & 1n) === 0n) {
      a >>= 1n;
      const bMod8 = ((b % 8n) + 8n) % 8n;
      if (bMod8 === 3n || bMod8 === 5n) {
        result = -result;
      }
    }
    [a, b] = [b, a];
    if (a % 4n === 3n && b % 4n === 3n) {
      result = -result;
    }
    a = a % b;
  }

  return b === 1n ? result : 0n;
}

/**
 * Cache of computed class numbers, keyed by the discriminant.
 *
 * Mirrors the `h_dict` local cache of Sage's
 * `discriminants_with_bounded_class_number`, but lives at module scope so that
 * repeated calls do not recompute anything.
 */
const h_dict: Map<bigint, bigint> = new Map();

/**
 * Compute the class number of the imaginary quadratic order of discriminant D
 * by counting reduced primitive positive-definite binary quadratic forms.
 *
 * A form (a, b, c) of discriminant `D = b^2 - 4ac < 0` with `a > 0` is reduced
 * iff `-a < b <= a <= c`, with `b >= 0` when `a == c`.  The class number is the
 * number of such forms that are primitive (`gcd(a, b, c) == 1`).
 *
 * @param D - a negative discriminant (D === 0 or 1 mod 4)
 * @returns the class number h(D)
 * @see Deviation: Sage calls PARI's `qfbclassno` (Shanks BSGS in the form class
 *   group); `qfbclassno` is not available in parigp-ts, so we count reduced
 *   forms instead.  The values are identical; the complexity is O(|D|) rather
 *   than O(|D|^(1/4)).
 */
function class_number(D: bigint): bigint {
  if (D >= 0n) {
    throw new ValueError(`D (= ${D}) must be negative`);
  }
  const cached = h_dict.get(D);
  if (cached !== undefined) {
    return cached;
  }

  const absD = -D;
  let count = 0n;
  // b runs over integers of the same parity as D with b^2 <= |D|/3
  let b = absD % 2n === 0n ? 0n : 1n;
  while (3n * b * b <= absD) {
    const q = (b * b + absD) / 4n; // = a*c
    let a = b === 0n ? 1n : b;
    if (a === 0n) a = 1n;
    while (a * a <= q) {
      if (q % a === 0n) {
        const c = q / a;
        if (_gcd3(a, b, c) === 1n) {
          // (a, b, c) is reduced; (a, -b, c) is a distinct reduced form
          // unless b === 0, a === b or a === c.
          count += b === 0n || a === b || a === c ? 1n : 2n;
        }
      }
      a += 1n;
    }
    b += 2n;
  }

  h_dict.set(D, count);
  return count;
}

/** gcd of three non-negative integers. */
function _gcd3(a: bigint, b: bigint, c: bigint): bigint {
  const g = (x: bigint, y: bigint): bigint => {
    let u = x < 0n ? -x : x;
    let v = y < 0n ? -y : y;
    while (v !== 0n) {
      [u, v] = [v, u % v];
    }
    return u;
  };
  return g(g(a, b), c);
}

/**
 * Compute the class number of a fundamental discriminant D.
 *
 * @param D - negative fundamental discriminant
 * @returns the class number h(D)
 */
function fundamental_class_number(D: bigint): bigint {
  return class_number(D);
}

/**
 * Return the class number h(f^2 * D0), given h(D0) = h0.
 *
 * @see Reference: sage/schemes/elliptic_curves/cm.py:OrderClassNumber
 */
function OrderClassNumber(D0: bigint, h0: bigint, f: bigint): bigint {
  if (f === 1n) return h0;

  // Check D0 is fundamental
  // (simplified check - full check would verify squarefree part)

  // Compute h(f^2 * D0) using the formula
  // h(D) = h(D0) * f * prod_{p|f} (1 - (D0/p)/p) / [O_{D0}^*: O_D^*]

  // Get prime factors of f
  const primes: bigint[] = [];
  let temp = f;
  let p = 2n;
  while (p * p <= temp) {
    if (temp % p === 0n) {
      primes.push(p);
      while (temp % p === 0n) {
        temp /= p;
      }
    }
    p++;
  }
  if (temp > 1n) {
    primes.push(temp);
  }

  // Compute the product
  let numerator = f;
  for (const prime of primes) {
    numerator /= prime;
    numerator *= prime - kronecker_symbol(D0, prime);
  }

  // Divide by the index of unit groups
  let units_index: bigint;
  if (D0 === -3n) {
    units_index = 3n;
  } else if (D0 === -4n) {
    units_index = 2n;
  } else {
    units_index = 1n;
  }

  return (h0 * numerator) / units_index;
}

/**
 * Return the Hilbert class polynomial for discriminant D.
 *
 * INPUT:
 * - D: negative integer congruent to 0 or 1 modulo 4
 * - algorithm: 'arb', 'sage', or 'magma' (default 'arb')
 *
 * OUTPUT: the Hilbert class polynomial H_D(x) in Z[x]
 *
 * The Hilbert class polynomial H_D is the minimal polynomial of the
 * j-invariant of the elliptic curve with CM by the order of discriminant D.
 *
 * EXAMPLES:
 * - hilbert_class_polynomial(-4) = x - 1728
 * - hilbert_class_polynomial(-7) = x + 3375
 * - hilbert_class_polynomial(-23) = x^3 + 3491750*x^2 - 5151296875*x + 12771880859375
 *
 * @see Reference: sage/schemes/elliptic_curves/cm.py:hilbert_class_polynomial
 */
export function hilbert_class_polynomial(
  D: bigint | number,
  algorithm?: 'arb' | 'sage' | 'magma'
): Polynomial<Integer> {
  const DVal = typeof D === 'number' ? BigInt(D) : D;

  // `cm.py` uses f-strings without a space: "D (=0) must be negative".
  if (DVal >= 0n) {
    throw new ValueError(`D (=${DVal}) must be negative`);
  }
  const Dmod4 = ((DVal % 4n) + 4n) % 4n;
  if (Dmod4 !== 0n && Dmod4 !== 1n) {
    throw new ValueError(`D (=${DVal}) must be a discriminant`);
  }
  if (algorithm === 'magma') {
    throw new NotImplementedError('Magma is not available');
  }

  // Delegate to PARI's `polclass` (`pari/src/basemath/polclass.c:1980-2108`),
  // ported in `parigp-ts/src/polmodular.ts` as `polclass0`.  This replaces a
  // nine-entry lookup table that raised NotImplementedError for every other
  // discriminant -- and whose hand-rolled `ZZ` object literal produced an
  // unprintable polynomial (`String(...)` gave `[object Object]*x + 3375`).
  const coeffs = polclass0(Number(DVal), INV_J, polmodular_db_init(INV_J));

  const R = new PolynomialRing<Integer>(ZZ, 'x');
  return R.__call__(coeffs.map((c) => new Integer(c)));
}

/**
 * Determine whether a polynomial is a Hilbert Class Polynomial.
 *
 * INPUT:
 * - f: a polynomial in Z[X]
 * - check_monic_irreducible: if True (default), check that f is monic and irreducible
 *
 * OUTPUT: D if f is the Hilbert Class Polynomial H_D, or 0 if not
 *
 * @see Reference: sage/schemes/elliptic_curves/cm.py:is_HCP
 */
export function is_HCP(_f: unknown, _check_monic_irreducible?: boolean): bigint {
  throw new NotImplementedError(
    'is_HCP requires polynomial root finding and elliptic curve endomorphism computation'
  );
}

/**
 * Return the CM j-invariants in a number field K.
 *
 * INPUT:
 * - K: a number field
 * - proof: if True (default), return proved results
 *
 * OUTPUT: a list of j-invariants in K that have complex multiplication
 *
 * @see Reference: sage/schemes/elliptic_curves/cm.py:cm_j_invariants
 */
export function cm_j_invariants(K: unknown, _proof?: boolean): unknown[] {
  // For K = QQ, return the known list
  if (
    K === 'QQ' ||
    (typeof K === 'object' &&
      K !== null &&
      'toString' in K &&
      (K as { toString: () => string }).toString() === 'Rational Field')
  ) {
    // Sage: sorted(j for D, f, j in cm_j_invariants_and_orders(K))
    return cm_j_invariants_QQ.map(([_, __, j]) => j).sort((a, b) => (a < b ? -1 : a > b ? 1 : 0));
  }

  throw new NotImplementedError(
    'cm_j_invariants for general number fields requires polynomial root finding'
  );
}

/**
 * Return pairs (j, O) of CM j-invariants and their CM orders.
 *
 * INPUT:
 * - K: a number field
 * - proof: if True (default), return proved results
 *
 * OUTPUT: a list of 3-tuples (D, f, j) where j is a CM j-invariant in K with
 * quadratic fundamental discriminant D and conductor f
 *
 * @see Reference: sage/schemes/elliptic_curves/cm.py:cm_j_invariants_and_orders
 */
export function cm_j_invariants_and_orders(
  K: unknown,
  _proof?: boolean
): Array<[bigint, bigint, bigint]> {
  // For K = QQ, return the known list as [D, f, j] tuples
  if (
    K === 'QQ' ||
    (typeof K === 'object' &&
      K !== null &&
      'toString' in K &&
      (K as { toString: () => string }).toString() === 'Rational Field')
  ) {
    // Sage returns 3-tuples (D, f, j)
    return cm_j_invariants_QQ.map(([d, f, j]) => [d, f, j] as [bigint, bigint, bigint]);
  }

  throw new NotImplementedError(
    'cm_j_invariants_and_orders for general number fields requires polynomial root finding'
  );
}

/**
 * Return all CM orders with a given class number h.
 *
 * INPUT:
 * - h: a positive integer (class number)
 * - proof: if True (default), return proved results
 *
 * OUTPUT: a list of orders in imaginary quadratic fields with class number h
 *
 * @see Reference: sage/schemes/elliptic_curves/cm.py:cm_orders
 */
export function cm_orders(h: bigint | number, _proof?: boolean): Array<[bigint, bigint]> {
  const hVal = typeof h === 'number' ? h : Number(h);

  if (hVal <= 0) {
    return [];
  }

  // Check if we have cached data
  if (hDf_dict.has(hVal)) {
    return hDf_dict.get(hVal)!;
  }

  // For h <= 100, compute using discriminants_with_bounded_class_number
  if (hVal <= 100) {
    const result = discriminants_with_bounded_class_number(BigInt(hVal));
    const orders = result.get(BigInt(hVal));
    return orders ?? [];
  }

  throw new NotImplementedError(`cm_orders for class number ${h} > 100 requires GRH bounds`);
}

/**
 * Return the largest fundamental discriminant with class number h.
 *
 * INPUT:
 * - h: a positive integer
 *
 * OUTPUT: a pair [|D|, n] where D is the largest (in absolute value)
 * fundamental discriminant with class number h and n is the number of such
 * discriminants; [0, 0] when h <= 0.
 *
 * @see Reference: sage/schemes/elliptic_curves/cm.py:largest_fundamental_disc_with_class_number
 */
export function largest_fundamental_disc_with_class_number(h: bigint | number): [bigint, bigint] {
  const hVal = typeof h === 'number' ? h : Number(h);

  if (hVal <= 0) {
    return [0n, 0n];
  }

  const entry = watkins_table.get(hVal);
  if (entry === undefined) {
    throw new NotImplementedError(
      `largest fundamental discriminant not available for class number ${h}`
    );
  }

  return [entry[0], BigInt(entry[1])];
}

/**
 * Return the largest discriminant with class number h.
 *
 * INPUT:
 * - h: a positive integer
 *
 * OUTPUT: a pair [|D|, n] where D is the largest (in absolute value)
 * discriminant with class number h and n is the number of such discriminants;
 * [0, 0] when h <= 0.
 *
 * @see Reference: sage/schemes/elliptic_curves/cm.py:largest_disc_with_class_number
 */
export function largest_disc_with_class_number(h: bigint | number): [bigint, bigint] {
  const hVal = typeof h === 'number' ? h : Number(h);

  if (hVal <= 0) {
    return [0n, 0n];
  }

  const entry = klaise_table.get(hVal);
  if (entry === undefined) {
    throw new NotImplementedError(`largest discriminant not available for class number ${h}`);
  }

  return [entry[0], BigInt(entry[1])];
}

/**
 * Return all discriminants with class number at most hmax.
 *
 * INPUT:
 * - hmax: maximum class number
 * - B: (optional) bound on discriminants
 * - proof: if True (default), return proved results
 *
 * OUTPUT: a dictionary mapping class numbers to lists of discriminants
 *
 * @see Reference: sage/schemes/elliptic_curves/cm.py:discriminants_with_bounded_class_number
 */
export function discriminants_with_bounded_class_number(
  hmax: bigint | number,
  B?: bigint | number,
  _proof?: boolean
): Map<bigint, Array<[bigint, bigint]>> {
  const hmaxVal = typeof hmax === 'number' ? BigInt(hmax) : hmax;
  const boundOption = B === undefined ? undefined : typeof B === 'number' ? BigInt(B) : B;

  // Easy case: everything up to hmax has already been computed and cached.
  // (cm.py:discriminants_with_bounded_class_number, "Easy case")
  const cachedMax = hDf_dict.size === 0 ? 0 : Math.max(...hDf_dict.keys());
  if (hDf_dict.size > 0 && hmaxVal <= BigInt(cachedMax)) {
    const T = new Map<bigint, Array<[bigint, bigint]>>();
    for (const [h, Dflist] of hDf_dict) {
      if (BigInt(h) > hmaxVal) continue;
      const list =
        boundOption === undefined
          ? [...Dflist]
          : Dflist.filter(([D0, f]) => (D0 < 0n ? -D0 : D0) * f * f <= boundOption);
      T.set(BigInt(h), list);
    }
    return T;
  }

  const result = new Map<bigint, Array<[bigint, bigint]>>();

  // Initialize result with empty arrays
  for (let h = 1n; h <= hmaxVal; h++) {
    result.set(h, []);
  }

  // Determine the bound B; `counts` is non-null exactly when B was *not*
  // supplied by the caller, in which case the result is complete and can both
  // be verified against Klaise's table and cached in hDf_dict.
  let bound: bigint;
  let counts: Map<bigint, number> | null = null;
  if (boundOption !== undefined) {
    bound = boundOption;
  } else if (hmaxVal <= 100n) {
    // Use the Klaise table to find the bound
    let maxDisc = 0n;
    counts = new Map();
    for (let h = 1n; h <= hmaxVal; h++) {
      const entry = klaise_table.get(Number(h));
      if (entry === undefined) {
        throw new NotImplementedError(`largest discriminant not available for class number ${h}`);
      }
      if (entry[0] > maxDisc) {
        maxDisc = entry[0];
      }
      counts.set(h, entry[1]);
    }
    bound = maxDisc;
  } else {
    throw new ValueError('if hmax > 100 you must specify a discriminant bound B');
  }

  // Iterate through discriminants
  for (let D = -3n; D >= -bound; D--) {
    // Check if D is a discriminant (D ≡ 0 or 1 mod 4)
    const Dmod4 = ((D % 4n) + 4n) % 4n;
    if (Dmod4 !== 0n && Dmod4 !== 1n) {
      continue;
    }

    // Find the fundamental part D0 and conductor f such that D = D0 * f^2
    let D0 = D;
    let f = 1n;

    // Find the squarefree part
    const absD = -D;
    let temp = absD;
    let product = 1n;

    for (let p = 2n; p * p <= temp; p++) {
      let count = 0n;
      while (temp % p === 0n) {
        temp /= p;
        count++;
      }
      if (count > 0n) {
        // If odd power, include in squarefree part
        if ((count & 1n) === 1n) {
          product *= p;
        }
        // Extract f
        f *= p ** (count / 2n);
      }
    }
    if (temp > 1n) {
      product *= temp;
    }

    // D0 is the fundamental discriminant
    D0 = -product;
    // Adjust D0 if needed (must be ≡ 0 or 1 mod 4)
    const D0mod4 = ((D0 % 4n) + 4n) % 4n;
    if (D0mod4 !== 0n && D0mod4 !== 1n) {
      D0 *= 4n;
      f /= 2n;
    }

    // Compute the class number: PARI's qfbclassno for fundamental
    // discriminants (cached in h_dict), OrderClassNumber for the rest.
    let h: bigint;
    if (f === 1n) {
      // D itself is fundamental
      h = fundamental_class_number(D0);
    } else {
      const cachedD = h_dict.get(D);
      if (cachedD !== undefined) {
        h = cachedD;
      } else {
        h = OrderClassNumber(D0, fundamental_class_number(D0), f);
        h_dict.set(D, h);
      }
    }

    // Add to result if class number is in range
    if (h >= 1n && h <= hmaxVal) {
      const list = result.get(h)!;
      list.push([D0, f]);
    }
  }

  // Sort each list by (|D0|, f)
  for (const [_, list] of result) {
    list.sort((a, b) => {
      const absD0_a = a[0] < 0n ? -a[0] : a[0];
      const absD0_b = b[0] < 0n ? -b[0] : b[0];
      if (absD0_a !== absD0_b) {
        return absD0_a < absD0_b ? -1 : 1;
      }
      return a[1] < b[1] ? -1 : a[1] > b[1] ? 1 : 0;
    });
  }

  if (counts !== null) {
    // 1. Check that we found the right number of discriminants
    for (const [h, list] of result) {
      const expected = counts.get(h);
      if (expected !== undefined && list.length !== expected) {
        throw new ValueError(
          `number of discriminants inconsistent with Klaise's table for h = ${h}: found ${list.length}, expected ${expected}`
        );
      }
    }
    // 2. Update the global cache
    for (const [h, list] of result) {
      hDf_dict.set(Number(h), list);
    }
  }

  return result;
}

/**
 * Determine whether j is a CM j-invariant.
 *
 * INPUT:
 * - j: an element of a number field
 * - algorithm: 'CremonaSutherland' (default) or other
 * - method: (deprecated) use algorithm instead
 *
 * OUTPUT: a pair [true, [D, f]] if j is the j-invariant of an order with
 * fundamental discriminant D and conductor f, otherwise [false, null]
 *
 * @see Reference: sage/schemes/elliptic_curves/cm.py:is_cm_j_invariant
 */
export function is_cm_j_invariant(
  j: unknown,
  _algorithm?: string,
  _method?: string
): [boolean, [bigint, bigint] | null] {
  // For rational integers, use the lookup table
  if (typeof j === 'bigint') {
    for (const [d, f, jInv] of cm_j_invariants_QQ) {
      if (jInv === j) {
        return [true, [d, f]];
      }
    }
    return [false, null];
  }

  if (typeof j === 'number') {
    return is_cm_j_invariant(BigInt(j), _algorithm, _method);
  }

  throw new NotImplementedError(
    'is_cm_j_invariant for general number field elements requires polynomial factorization'
  );
}
