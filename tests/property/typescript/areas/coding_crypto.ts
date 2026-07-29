/**
 * sagemath-ts side of the `coding_crypto` property-test area.
 *
 * Cases: tests/property/cases/coding_crypto.cases.json
 * SageMath counterpart: tests/property/python/areas/coding_crypto.py
 *
 * Covers `sage.coding.*` (Reed-Muller, BCH, generalized Reed-Solomon, Goppa)
 * and `sage.crypto.*` (S-boxes, Boolean functions, LWE parameter derivation,
 * hard-lattice generation).
 *
 * Conventions shared with the SageMath side (both files must agree exactly):
 *
 * - Every function returns an **already formatted string**, so the generic
 *   `formatResult` / `format_result` in the runners never sees a shape it has
 *   to guess about (its `Factorization` heuristic would otherwise mangle a
 *   2 x 2 difference distribution table).
 * - Lists are rendered `[a, b, c]`, matrices `[[a, b], [c, d]]`, booleans
 *   `True` / `False` -- i.e. Python's `str()`.
 * - A function that is *supposed* to raise on both sides is wrapped in
 *   `guard()`, which renders any exception as the single token `ERROR`.
 *   `compare.ts` already scores "both runners raised" as a pass, so this only
 *   adds strength (one side raising while the other returns a value now fails);
 *   it deliberately does not compare exception messages.
 * - Finite field elements of GF(p^k) are passed in and out as **integers**:
 *   the integer whose base-p digits are the coefficients of the element in the
 *   polynomial basis (SageMath's `F.from_integer` / `x.to_integer()`, the port's
 *   `F.fromInteger` / `x.integer_representation()`).  That keeps the case JSON
 *   free of field-specific syntax and pins the *representation* too.
 *
 * Port defects this area currently pins (28 of 926 cases fail; every expected
 * value came from running `sage`, so the fix belongs in the port, never in the
 * case):
 *
 * | symptom | cases | root cause |
 * |---|---|---|
 * | `sbox_*` on `SBox([0,0,0,0])`, `SBox([0, 2**49])` | 15 | `crypto/sbox.ts:228` uses `max(1, ceil(log2(max+1)))` for the output size; SageMath uses `ZZ(max(S)).nbits()` (`sbox.pyx:208`) -- 0 for the constant S-box, and exact above 2^53 |
 * | `sbox_is_involution` on non-permutations | 5 | `crypto/sbox.ts:383` returns `false`; SageMath's is `self == self.inverse()` (`sbox.pyx:1914`) and `inverse` raises `TypeError` (`:1811`) |
 * | `rm_parameters_formula(2, 60)`, `(0, 63)` | 2 | `coding/reed_muller_code.ts:232` keeps the length in a JS `number` |
 * | `goppa_generator_matrix` | 4 | `coding/goppa_code.ts:451` returns a non-echelon basis; SageMath returns `from_parity_check_matrix(H).generator_matrix()` (`goppa_code.py:434-437`) |
 * | `lattice_gen_random` with q >= 2^31 | 2 | `crypto/lattice.ts:196` always draws one 31-bit `c_random()`; SageMath only does that for the small-modulus dense matrix templates |
 */

import { BCHCode } from '../../../../packages/sagemath-ts/src/coding/bch_code.js';
import { GoppaCode } from '../../../../packages/sagemath-ts/src/coding/goppa_code.js';
import { ReedMullerCode } from '../../../../packages/sagemath-ts/src/coding/reed_muller_code.js';
import { ReedSolomonCode } from '../../../../packages/sagemath-ts/src/coding/reed_solomon.js';
import { BooleanFunction } from '../../../../packages/sagemath-ts/src/crypto/boolean_function.js';
import { gen_lattice } from '../../../../packages/sagemath-ts/src/crypto/lattice.js';
import {
  LindnerPeikert,
  Regev,
  RingLindnerPeikert,
} from '../../../../packages/sagemath-ts/src/crypto/lwe.js';
import {
  SBox,
  feistel_construction,
  misty_construction,
} from '../../../../packages/sagemath-ts/src/crypto/sbox.js';
import { GF } from '../../../../packages/sagemath-ts/src/rings/finite_rings/finite_field_constructor.js';
import { GFExtended } from '../../../../packages/sagemath-ts/src/rings/finite_rings/finite_field_extension.js';
import { Polynomial } from '../../../../packages/sagemath-ts/src/rings/polynomial/polynomial_element.js';
import { PolynomialRing } from '../../../../packages/sagemath-ts/src/rings/polynomial/polynomial_ring.js';

// ---------------------------------------------------------------------------
// formatting helpers (mirrored in the Python area module)
// ---------------------------------------------------------------------------

/** Render a flat list the way Python's `str(list)` does. */
function lst(xs: Iterable<unknown>): string {
  return `[${[...xs].map((x) => String(x)).join(', ')}]`;
}

/** Render a list of lists the way Python's `str(list)` does. */
function mat(rows: Iterable<Iterable<unknown>>): string {
  return `[${[...rows].map((r) => lst(r)).join(', ')}]`;
}

/** Render a boolean the way Python does. */
function bool_(b: boolean): string {
  return b ? 'True' : 'False';
}

/** Run `f`, collapsing any exception to the single token `ERROR`. */
function guard(f: () => string): string {
  try {
    return f();
  } catch {
    return 'ERROR';
  }
}

const num = (x: bigint | number): number => Number(x);
const nums = (xs: (bigint | number)[]): number[] => xs.map(num);

// ---------------------------------------------------------------------------
// sage.crypto.sbox
// ---------------------------------------------------------------------------

function sbox(S: bigint[], bigEndian: bigint = 1n): SBox {
  return new SBox(nums(S), { big_endian: num(bigEndian) !== 0 });
}

function matrixRows<T extends { toString(): string }>(M: {
  nrows: number;
  ncols: number;
  get(i: number, j: number): T;
}): string[][] {
  const rows: string[][] = [];
  for (let i = 0; i < M.nrows; i++) {
    const row: string[] = [];
    for (let j = 0; j < M.ncols; j++) {
      row.push(M.get(i, j).toString());
    }
    rows.push(row);
  }
  return rows;
}

// ---------------------------------------------------------------------------
// sage.crypto.boolean_function
// ---------------------------------------------------------------------------

function bf(T: bigint[]): BooleanFunction {
  return new BooleanFunction(nums(T) as (0 | 1)[]);
}

// ---------------------------------------------------------------------------
// sage.coding helpers
// ---------------------------------------------------------------------------

/** GF(p) or GF(p^k); `fromInt` maps SageMath's integer representation in. */
// biome-ignore lint/suspicious/noExplicitAny: the two field classes are structurally different
function field(q: bigint): any {
  return GFExtended(q);
}

// biome-ignore lint/suspicious/noExplicitAny: element type depends on the field
function fromInt(F: any, i: bigint): any {
  if (typeof F.fromInteger === 'function') {
    return F.fromInteger(i);
  }
  return F.__call__(i);
}

// biome-ignore lint/suspicious/noExplicitAny: element type depends on the field
function toInt(F: any, x: any): bigint {
  if (typeof x.integer_representation === 'function') {
    return x.integer_representation();
  }
  return x.toBigInt();
}

/** Codeword of the cyclic code generated by `g`: the coefficients of m(x)g(x). */
// biome-ignore lint/suspicious/noExplicitAny: polynomial element type depends on the field
function cyclicCodeword(F: any, g: any, msg: bigint[], n: number): any[] {
  const R = g.parent as PolynomialRing<never>;
  const m = new Polynomial(
    msg.map((c) => fromInt(F, c)),
    R as never
  );
  const prod = m.mul(g);
  const out = [];
  for (let i = 0; i < n; i++) {
    out.push(prod.getCoeff(i));
  }
  return out;
}

// ---------------------------------------------------------------------------

export const functions = {
  // -- sage.crypto.sbox ----------------------------------------------------

  /** `(input_size(), output_size())` -- SageMath uses `ZZ(max(S)).nbits()`. */
  sbox_sizes: (S: bigint[]) => guard(() => lst([sbox(S).input_size(), sbox(S).output_size()])),

  sbox_ddt: (S: bigint[]) => guard(() => mat(matrixRows(sbox(S).difference_distribution_table()))),

  sbox_differential_uniformity: (S: bigint[]) =>
    guard(() => String(sbox(S).differential_uniformity())),

  sbox_is_apn: (S: bigint[]) => guard(() => bool_(sbox(S).is_apn())),

  /** `scale` index: 0 bias, 1 correlation, 2 absolute_bias (default), 3 fourier_coefficient. */
  sbox_lat: (S: bigint[], scale: bigint) =>
    guard(() => {
      const scales = ['bias', 'correlation', 'absolute_bias', 'fourier_coefficient'] as const;
      const sc = scales[num(scale)];
      if (sc === undefined) {
        throw new Error('bad scale index');
      }
      return mat(matrixRows(sbox(S).linear_approximation_table(sc)));
    }),

  sbox_linearity: (S: bigint[]) => guard(() => String(sbox(S).linearity())),
  sbox_max_linear_bias: (S: bigint[]) => guard(() => String(sbox(S).max_linear_bias())),
  sbox_nonlinearity: (S: bigint[]) => guard(() => String(sbox(S).nonlinearity())),
  sbox_max_degree: (S: bigint[]) => guard(() => String(sbox(S).max_degree())),
  sbox_min_degree: (S: bigint[]) => guard(() => String(sbox(S).min_degree())),
  sbox_is_permutation: (S: bigint[]) => guard(() => bool_(sbox(S).is_permutation())),
  sbox_is_involution: (S: bigint[]) => guard(() => bool_(sbox(S).is_involution())),
  sbox_is_balanced: (S: bigint[]) => guard(() => bool_(sbox(S).is_balanced())),
  sbox_inverse: (S: bigint[]) => guard(() => lst(sbox(S).inverse().toArray())),
  sbox_derivative: (S: bigint[], u: bigint) =>
    guard(() => lst(sbox(S).derivative(num(u)).toArray())),
  sbox_component_function: (S: bigint[], b: bigint) =>
    guard(() => lst(sbox(S).component_function_table(num(b)))),
  sbox_fixed_points: (S: bigint[]) => guard(() => lst(sbox(S).fixed_points())),
  sbox_differential_branch_number: (S: bigint[]) =>
    guard(() => String(sbox(S).differential_branch_number())),
  sbox_linear_branch_number: (S: bigint[]) => guard(() => String(sbox(S).linear_branch_number())),

  /** `S(x)` on bit vectors, exercising the `big_endian` flag. */
  sbox_call_bits: (S: bigint[], x: bigint, bigEndian: bigint) =>
    guard(() => {
      const s = sbox(S, bigEndian);
      return lst(s.call(s.to_bits(num(x), s.input_size())) as number[]);
    }),

  sbox_to_bits: (S: bigint[], x: bigint, n: bigint, bigEndian: bigint) =>
    guard(() => lst(sbox(S, bigEndian).to_bits(num(x), num(n)))),

  sbox_from_bits: (S: bigint[], bits: bigint[], bigEndian: bigint) =>
    guard(() => String(sbox(S, bigEndian).from_bits(nums(bits)))),

  /** Feistel construction from 2 or 3 round S-boxes (`sbox.pyx:1990`). */
  sbox_feistel: (a: bigint[], b: bigint[], c: bigint[]) =>
    guard(() => {
      const boxes = [sbox(a), sbox(b)];
      if (c.length > 0) {
        boxes.push(sbox(c));
      }
      return lst(feistel_construction(boxes).toArray());
    }),

  /** MISTY construction (`sbox.pyx:2010`) -- H108 swapped the halves. */
  sbox_misty: (a: bigint[], b: bigint[], c: bigint[]) =>
    guard(() => {
      const boxes = [sbox(a), sbox(b)];
      if (c.length > 0) {
        boxes.push(sbox(c));
      }
      return lst(misty_construction(boxes).toArray());
    }),

  sbox_misty_stats: (a: bigint[], b: bigint[], c: bigint[]) =>
    guard(() => {
      const boxes = [sbox(a), sbox(b)];
      if (c.length > 0) {
        boxes.push(sbox(c));
      }
      const S = misty_construction(boxes);
      return lst([S.differential_uniformity(), S.linearity()]);
    }),

  // -- sage.crypto.boolean_function ---------------------------------------

  bf_nvariables: (T: bigint[]) => guard(() => String(bf(T).nvariables())),
  bf_truth_table_hex: (T: bigint[]) => guard(() => String(bf(T).truthTable('hex'))),

  /** Build from the hex string `value` zero-padded to `width` characters. */
  bf_from_hex: (value: bigint, width: bigint) =>
    guard(() => {
      const hex = value.toString(16).padStart(num(width), '0');
      return lst(new BooleanFunction(hex).truthTable('int') as number[]);
    }),

  bf_walsh: (T: bigint[]) => guard(() => lst(bf(T).walshHadamardTransform())),
  bf_absolute_walsh_spectrum: (T: bigint[]) =>
    guard(() => {
      const spec = [...bf(T).absoluteWalshSpectrum().entries()].sort((x, y) => x[0] - y[0]);
      return lst(spec.map(([k, v]) => `(${k}, ${v})`));
    }),
  bf_nonlinearity: (T: bigint[]) => guard(() => String(bf(T).nonlinearity())),
  bf_is_bent: (T: bigint[]) => guard(() => bool_(bf(T).isBent())),
  bf_is_balanced: (T: bigint[]) => guard(() => bool_(bf(T).isBalanced())),
  bf_correlation_immunity: (T: bigint[]) => guard(() => String(bf(T).correlationImmunity())),
  bf_resiliency_order: (T: bigint[]) => guard(() => String(bf(T).resiliencyOrder())),
  bf_autocorrelation: (T: bigint[]) => guard(() => lst(bf(T).autocorrelation())),
  bf_absolute_indicator: (T: bigint[]) => guard(() => String(bf(T).absoluteIndicator())),
  bf_sum_of_square_indicator: (T: bigint[]) => guard(() => String(bf(T).sumOfSquareIndicator())),
  bf_is_plateaued: (T: bigint[]) => guard(() => bool_(bf(T).isPlateaued())),
  bf_algebraic_degree: (T: bigint[]) => guard(() => String(bf(T).algebraicDegree())),
  bf_algebraic_immunity: (T: bigint[]) => guard(() => String(bf(T).algebraicImmunity())),
  bf_anf_coefficients: (T: bigint[]) => guard(() => lst(bf(T).algebraicNormalFormCoefficients())),
  bf_derivative: (T: bigint[], u: bigint) =>
    guard(() => lst(bf(T).derivative(num(u)).truthTable('int') as number[])),
  bf_is_linear_structure: (T: bigint[], u: bigint) =>
    guard(() => bool_(bf(T).isLinearStructure(num(u)))),
  bf_is_linear_structure_vec: (T: bigint[], u: bigint[]) =>
    guard(() => bool_(bf(T).isLinearStructure(nums(u)))),
  bf_has_linear_structure: (T: bigint[]) => guard(() => bool_(bf(T).hasLinearStructure())),
  bf_is_symmetric: (T: bigint[]) => guard(() => bool_(bf(T).isSymmetric())),
  bf_complement: (T: bigint[]) =>
    guard(() => lst(bf(T).complement().truthTable('int') as number[])),
  bf_add: (A: bigint[], B: bigint[]) =>
    guard(() => lst(bf(A).add(bf(B)).truthTable('int') as number[])),
  bf_mul: (A: bigint[], B: bigint[]) =>
    guard(() => lst(bf(A).mul(bf(B)).truthTable('int') as number[])),
  bf_concatenate: (A: bigint[], B: bigint[]) =>
    guard(() => lst(bf(A).concatenate(bf(B)).truthTable('int') as number[])),
  /** The `BooleanFunction(n)` constructor: the zero function on n variables. */
  bf_zero_function: (n: bigint) =>
    guard(() => {
      const f = new BooleanFunction(num(n));
      return lst([f.nvariables(), ...(f.truthTable('int') as number[])]);
    }),

  // -- sage.coding.reed_muller_code ---------------------------------------

  rm_parameters: (r: bigint, m: bigint) =>
    guard(() => {
      const C = new ReedMullerCode(r, m);
      return lst([C.length(), C.dimension(), C.minimum_distance()]);
    }),

  /**
   * Same three numbers as `rm_parameters`, for lengths SageMath 10.3 cannot
   * instantiate (it segfaults building the ambient GF(2)^(2^m)); the SageMath
   * side then evaluates the formulas upstream's own class uses.  M124: the port
   * used to compute the length with a 32-bit `1 << m`.
   */
  rm_parameters_formula: (r: bigint, m: bigint) =>
    guard(() => {
      const C = new ReedMullerCode(r, m);
      return lst([C.length(), C.dimension(), C.minimum_distance()]);
    }),

  rm_generator_matrix: (r: bigint, m: bigint) =>
    guard(() => mat(matrixRows(new ReedMullerCode(r, m).generator_matrix()))),

  rm_encode: (r: bigint, m: bigint, msg: bigint[]) =>
    guard(() => lst(new ReedMullerCode(r, m).encode(nums(msg)).map((x) => x.value))),

  /**
   * Encode, flip the given positions, decode back to the message.
   *
   * Note the naming: `ReedMullerCode.decode` (and `ReedSolomonCode.decode`)
   * return the *message*, i.e. they are SageMath's `decode_to_message`, not its
   * `decode_to_code`; re-encoding the result gives the corrected codeword.
   */
  rm_decode_to_message: (r: bigint, m: bigint, msg: bigint[], errs: bigint[]) =>
    guard(() => {
      const C = new ReedMullerCode(r, m);
      const cw = C.encode(nums(msg)).map((x) => Number(x.value));
      for (const e of nums(errs)) {
        cw[e] = 1 - cw[e]!;
      }
      return lst(C.decode(cw).map((x) => x.value));
    }),

  /** Same, but re-encoded: the corrected *codeword*, as SageMath's `decode_to_code`. */
  rm_decode_to_code: (r: bigint, m: bigint, msg: bigint[], errs: bigint[]) =>
    guard(() => {
      const C = new ReedMullerCode(r, m);
      const cw = C.encode(nums(msg)).map((x) => Number(x.value));
      for (const e of nums(errs)) {
        cw[e] = 1 - cw[e]!;
      }
      return lst(C.encode(C.decode(cw)).map((x) => x.value));
    }),

  // -- sage.coding.bch_code -----------------------------------------------

  bch_generator_polynomial: (q: bigint, n: bigint, delta: bigint, b: bigint, l: bigint) =>
    guard(() => {
      const F = field(q);
      const C = new BCHCode(n, delta, F, b, l);
      const g = C.generator_polynomial();
      const out: string[] = [];
      for (let i = 0; i <= g.degree(); i++) {
        out.push(g.getCoeff(i).toString());
      }
      return lst(out);
    }),

  bch_defining_set: (q: bigint, n: bigint, delta: bigint, b: bigint, l: bigint) =>
    guard(() =>
      lst([...new BCHCode(n, delta, field(q), b, l).defining_set()].sort((x, y) => x - y))
    ),

  bch_dimension: (q: bigint, n: bigint, delta: bigint, b: bigint, l: bigint) =>
    guard(() => String(new BCHCode(n, delta, field(q), b, l).dimension())),

  /**
   * m(x)*g(x), corrupted at `errs` (position, value pairs flattened), decoded.
   * Prime base fields only, so the symbols are plain integers.
   */
  bch_decode: (
    q: bigint,
    n: bigint,
    delta: bigint,
    b: bigint,
    l: bigint,
    msg: bigint[],
    errs: bigint[]
  ) =>
    guard(() => {
      const F = field(q);
      const C = new BCHCode(n, delta, F, b, l);
      const g = C.generator_polynomial();
      const cw = cyclicCodeword(F, g, msg, num(n));
      const v = [...cw];
      for (let i = 0; i < errs.length; i += 2) {
        const pos = num(errs[i]!);
        v[pos] = v[pos].add(fromInt(F, errs[i + 1]!));
      }
      return lst(C.decode(v).map((x) => toInt(F, x)));
    }),

  // -- sage.coding.grs_code (Reed-Solomon) --------------------------------

  /** GRS code over GF(p) with evaluation points 0, 1, ..., n-1 and unit multipliers. */
  grs_encode: (p: bigint, n: bigint, k: bigint, msg: bigint[]) =>
    guard(() => {
      const F = GF(p);
      const pts = Array.from({ length: num(n) }, (_, i) => F.__call__(BigInt(i)));
      const C = new ReedSolomonCode(F, n, k, pts);
      return lst(C.encode(msg.map((c) => F.__call__(c))).map((x) => x.toBigInt()));
    }),

  grs_parity_column_multipliers: (p: bigint, n: bigint, k: bigint) =>
    guard(() => {
      const F = GF(p);
      const pts = Array.from({ length: num(n) }, (_, i) => F.__call__(BigInt(i)));
      const C = new ReedSolomonCode(F, n, k, pts);
      return lst(C.parity_column_multipliers().map((x) => x.toBigInt()));
    }),

  grs_syndrome: (p: bigint, n: bigint, k: bigint, word: bigint[]) =>
    guard(() => {
      const F = GF(p);
      const pts = Array.from({ length: num(n) }, (_, i) => F.__call__(BigInt(i)));
      const C = new ReedSolomonCode(F, n, k, pts);
      return lst(C.syndrome(word.map((c) => F.__call__(c))).map((x) => x.toBigInt()));
    }),

  grs_decode_to_message: (p: bigint, n: bigint, k: bigint, msg: bigint[], errs: bigint[]) =>
    guard(() => {
      const F = GF(p);
      const pts = Array.from({ length: num(n) }, (_, i) => F.__call__(BigInt(i)));
      const C = new ReedSolomonCode(F, n, k, pts);
      const cw = C.encode(msg.map((c) => F.__call__(c)));
      for (let i = 0; i < errs.length; i += 2) {
        const pos = num(errs[i]!);
        cw[pos] = cw[pos]!.add(F.__call__(errs[i + 1]!));
      }
      return lst(C.decode(cw).map((x) => x.toBigInt()));
    }),

  grs_decode_to_code: (p: bigint, n: bigint, k: bigint, msg: bigint[], errs: bigint[]) =>
    guard(() => {
      const F = GF(p);
      const pts = Array.from({ length: num(n) }, (_, i) => F.__call__(BigInt(i)));
      const C = new ReedSolomonCode(F, n, k, pts);
      const cw = C.encode(msg.map((c) => F.__call__(c)));
      for (let i = 0; i < errs.length; i += 2) {
        const pos = num(errs[i]!);
        cw[pos] = cw[pos]!.add(F.__call__(errs[i + 1]!));
      }
      return lst(C.encode(C.decode(cw)).map((x) => x.toBigInt()));
    }),

  // -- sage.coding.goppa_code ---------------------------------------------

  goppa_parity_check_matrix: (q: bigint, gcoeffs: bigint[], support: bigint[]) =>
    guard(() => {
      const F = field(q);
      const R = new PolynomialRing(F, 'x');
      const g = new Polynomial(
        gcoeffs.map((c) => fromInt(F, c)),
        R
      );
      const L = support.map((i) => fromInt(F, i));
      const C = new GoppaCode(g, L);
      return mat(matrixRows(C.parity_check_matrix()));
    }),

  goppa_dimension: (q: bigint, gcoeffs: bigint[], support: bigint[]) =>
    guard(() => {
      const F = field(q);
      const R = new PolynomialRing(F, 'x');
      const g = new Polynomial(
        gcoeffs.map((c) => fromInt(F, c)),
        R
      );
      const C = new GoppaCode(
        g,
        support.map((i) => fromInt(F, i))
      );
      return lst([C.length(), C.dimension()]);
    }),

  goppa_distance_bound: (q: bigint, gcoeffs: bigint[], support: bigint[]) =>
    guard(() => {
      const F = field(q);
      const R = new PolynomialRing(F, 'x');
      const g = new Polynomial(
        gcoeffs.map((c) => fromInt(F, c)),
        R
      );
      return String(
        new GoppaCode(
          g,
          support.map((i) => fromInt(F, i))
        ).distance_bound()
      );
    }),

  goppa_generator_matrix: (q: bigint, gcoeffs: bigint[], support: bigint[]) =>
    guard(() => {
      const F = field(q);
      const R = new PolynomialRing(F, 'x');
      const g = new Polynomial(
        gcoeffs.map((c) => fromInt(F, c)),
        R
      );
      const C = new GoppaCode(
        g,
        support.map((i) => fromInt(F, i))
      );
      return mat(matrixRows(C.generator_matrix()));
    }),

  /** `codeword` is a genuine codeword (checked by both sides via H); flip `errs` and decode. */
  goppa_decode: (
    q: bigint,
    gcoeffs: bigint[],
    support: bigint[],
    codeword: bigint[],
    errs: bigint[]
  ) =>
    guard(() => {
      const F = field(q);
      const R = new PolynomialRing(F, 'x');
      const g = new Polynomial(
        gcoeffs.map((c) => fromInt(F, c)),
        R
      );
      const C = new GoppaCode(
        g,
        support.map((i) => fromInt(F, i))
      );
      const word = nums(codeword);
      for (const e of nums(errs)) {
        word[e] = 1 - word[e]!;
      }
      const decoded = C.decode(word.map((b) => fromInt(F, BigInt(b))));
      return lst(decoded.map((x) => toInt(F, x)));
    }),

  // -- sage.crypto.lwe -----------------------------------------------------

  /** `(q, sigma, m, secret_dist)` of Regev's parameter set. */
  lwe_regev_params: (n: bigint) =>
    guard(() => {
      const l = new Regev(n);
      return lst([l.K.order, l.D.sigma.toFixed(6), l.m === null ? 'None' : l.m, l.secret_dist]);
    }),

  lwe_lindner_peikert_params: (n: bigint) =>
    guard(() => {
      const l = new LindnerPeikert(n);
      return lst([l.K.order, l.D.sigma.toFixed(6), l.m === null ? 'None' : l.m, l.secret_dist]);
    }),

  lwe_ring_lindner_peikert_params: (n: bigint) =>
    guard(() => {
      const l = new RingLindnerPeikert(n);
      return lst([l.K.order, l.D.sigma.toFixed(6), l.m === null ? 'None' : l.m]);
    }),

  // -- sage.crypto.lattice -------------------------------------------------

  lattice_gen_modular: (n: bigint, m: bigint, q: bigint, seed: bigint) =>
    guard(() => mat(gen_lattice({ type: 'modular', n, m, q, seed: num(seed) }) as bigint[][])),

  lattice_gen_modular_dual: (n: bigint, m: bigint, q: bigint, seed: bigint) =>
    guard(() =>
      mat(gen_lattice({ type: 'modular', n, m, q, seed: num(seed), dual: true }) as bigint[][])
    ),

  lattice_gen_random: (m: bigint, q: bigint, seed: bigint) =>
    guard(() => mat(gen_lattice({ type: 'random', n: 1n, m, q, seed: num(seed) }) as bigint[][])),

  /**
   * Shape and determinant of a generated basis: `det = q^n` (primal) and
   * `q^(m-n)` (dual).  H69/H107 swapped the two blocks, which this pins
   * without depending on the random entries.
   */
  lattice_gen_invariants: (n: bigint, m: bigint, q: bigint, seed: bigint, dual: bigint) =>
    guard(() => {
      const B = gen_lattice({
        type: 'modular',
        n,
        m,
        q,
        seed: num(seed),
        dual: num(dual) !== 0,
      }) as bigint[][];
      let det = detBigint(B);
      if (det < 0n) {
        det = -det;
      }
      let qRows = 0;
      for (const row of B) {
        if (row.every((x) => x === 0n || x === q || x === -q)) {
          qRows++;
        }
      }
      return lst([B.length, B[0]!.length, det, qRows]);
    }),

  lattice_gen_cyclotomic_invariants: (n: bigint, m: bigint, q: bigint, seed: bigint) =>
    guard(() => {
      const B = gen_lattice({ type: 'cyclotomic', n, m, q, seed: num(seed) }) as bigint[][];
      let det = detBigint(B);
      if (det < 0n) {
        det = -det;
      }
      let qRows = 0;
      for (const row of B) {
        if (row.every((x) => x === 0n || x === q || x === -q)) {
          qRows++;
        }
      }
      return lst([B.length, B[0]!.length, det, qRows]);
    }),
};

/** Exact integer determinant by fraction-free (Bareiss) elimination. */
function detBigint(M: bigint[][]): bigint {
  const n = M.length;
  const A = M.map((r) => [...r]);
  let sign = 1n;
  let prev = 1n;
  for (let k = 0; k < n - 1; k++) {
    if (A[k]![k] === 0n) {
      let swap = -1;
      for (let i = k + 1; i < n; i++) {
        if (A[i]![k] !== 0n) {
          swap = i;
          break;
        }
      }
      if (swap === -1) {
        return 0n;
      }
      const tmp = A[k]!;
      A[k] = A[swap]!;
      A[swap] = tmp;
      sign = -sign;
    }
    for (let i = k + 1; i < n; i++) {
      for (let j = k + 1; j < n; j++) {
        A[i]![j] = (A[i]![j]! * A[k]![k]! - A[i]![k]! * A[k]![j]!) / prev;
      }
      A[i]![k] = 0n;
    }
    prev = A[k]![k]!;
  }
  return sign * A[n - 1]![n - 1]!;
}
