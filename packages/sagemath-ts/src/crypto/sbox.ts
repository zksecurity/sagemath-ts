/**
 * @module sage/crypto/sbox
 * @description S-Box cryptographic analysis
 *
 * Port of: sage/crypto/sbox.pyx
 *
 * Implements substitution boxes (S-boxes) which are fundamental components
 * of symmetric key cryptography. Provides tools for computing cryptographic
 * properties like Difference Distribution Tables (DDT), Linear Approximation
 * Tables (LAT), differential uniformity, and more.
 *
 * @example
 * ```typescript
 * import { SBox } from 'sagemath-ts/crypto/sbox';
 *
 * // Create the PRESENT S-box
 * const S = new SBox([12, 5, 6, 11, 9, 0, 10, 13, 3, 14, 15, 8, 4, 7, 1, 2]);
 *
 * // Evaluate S-box
 * console.log(S.call(1)); // 5
 *
 * // Compute DDT
 * const ddt = S.difference_distribution_table();
 * console.log(S.differential_uniformity()); // Maximum non-trivial DDT entry
 *
 * // Compute LAT
 * const lat = S.linear_approximation_table();
 * console.log(S.max_linear_bias()); // Maximum absolute LAT entry
 * ```
 */

import { TypeError as SageTypeError, ValueError } from '../errors.js';
import { Matrix } from '../matrix/matrix_generic.js';
import { Integer, ZZ } from '../rings/integer_ring.js';
import { BooleanFunction } from './boolean_function.js';

/**
 * Helper to compute Hamming weight (population count) of a number.
 */
function hammingWeight(x: number): number {
  let count = 0;
  while (x > 0) {
    count += x & 1;
    x >>>= 1;
  }
  return count;
}

/**
 * Helper to compute dot product in GF(2) (XOR of ANDed bits).
 * Returns 1 if odd number of 1s in (a & b), else 0.
 */
function gf2DotProduct(a: number, b: number): number {
  return hammingWeight(a & b) & 1;
}

/**
 * Interface for ring elements used in the Matrix class.
 * We create a simple integer ring for use with our Matrix.
 */
class IntegerElement {
  constructor(public readonly value: bigint) {}

  add(other: IntegerElement): IntegerElement {
    return new IntegerElement(this.value + other.value);
  }

  sub(other: IntegerElement): IntegerElement {
    return new IntegerElement(this.value - other.value);
  }

  mul(other: IntegerElement): IntegerElement {
    return new IntegerElement(this.value * other.value);
  }

  neg(): IntegerElement {
    return new IntegerElement(-this.value);
  }

  isZero(): boolean {
    return this.value === 0n;
  }

  eq(other: IntegerElement): boolean {
    return this.value === other.value;
  }

  toString(): string {
    return this.value.toString();
  }
}

/**
 * Simple integer ring for matrix operations.
 */
const IntegerRingForMatrix = {
  zero: () => new IntegerElement(0n),
  one: () => new IntegerElement(1n),
  __call__: (x: bigint | number) => new IntegerElement(BigInt(x)),
};

/**
 * A substitution box (S-box) for cryptographic analysis.
 *
 * An S-box is one of the basic components of symmetric key cryptography.
 * In general, an S-box takes m input bits and transforms them into n output bits.
 * This is called an m x n S-box and is often implemented as a lookup table.
 *
 * This class allows algebraic treatment and determination of various
 * cryptographic properties including:
 * - Difference Distribution Table (DDT)
 * - Linear Approximation Table (LAT)
 * - Differential uniformity
 * - Linearity and nonlinearity
 * - APN (Almost Perfect Nonlinear) property
 *
 * @example
 * ```typescript
 * // Create the PRESENT S-box
 * const S = new SBox([12, 5, 6, 11, 9, 0, 10, 13, 3, 14, 15, 8, 4, 7, 1, 2]);
 *
 * // Evaluate
 * console.log(S.call(1)); // 5
 * console.log(S.call([0, 0, 0, 1])); // [0, 1, 0, 1]
 *
 * // Check properties
 * console.log(S.is_permutation()); // true
 * console.log(S.differential_uniformity()); // 4
 * ```
 */
export class SBox {
  private readonly _S_list: number[];
  private readonly _m: number; // input size in bits
  private readonly _n: number; // output size in bits
  private readonly _big_endian: boolean;

  // Cached tables
  private _ddt_cache: Matrix<IntegerElement> | null = null;
  private _lat_cache: Matrix<IntegerElement> | null = null;

  /**
   * Construct a substitution box (S-box) for a given lookup table.
   *
   * @param S - The lookup table defining the S-box. Can be an array of integers
   *            where S[i] is the output for input i.
   * @param options - Configuration options
   * @param options.big_endian - If true, bits are ordered big-endian (MSB first).
   *                             Default is true, consistent with cryptographic literature.
   *
   * @throws {SageTypeError} If no lookup table is provided
   * @throws {SageTypeError} If lookup table length is not a power of 2
   *
   * @example
   * ```typescript
   * // 3-bit S-box where (0,0,1) maps to (1,1,1)
   * const S = new SBox([7, 6, 0, 4, 2, 5, 1, 3]);
   * console.log(S.call(0)); // 7
   *
   * // Little-endian S-box
   * const S_le = new SBox([7, 6, 0, 4, 2, 5, 1, 3], { big_endian: false });
   * ```
   */
  constructor(S: number[], options?: { big_endian?: boolean }) {
    if (!S || S.length === 0) {
      throw new SageTypeError('no lookup table provided');
    }

    // Check that length is a power of 2
    const len = S.length;
    if ((len & (len - 1)) !== 0) {
      throw new SageTypeError('lookup table length is not a power of 2');
    }

    this._S_list = [...S];
    this._m = Math.log2(len);
    this._n = Math.max(1, Math.ceil(Math.log2(Math.max(...S) + 1)));
    this._big_endian = options?.big_endian ?? true;
  }

  /**
   * Return the length of input bit strings.
   */
  get length(): number {
    return this._m;
  }

  /**
   * Return the input size of this S-Box (in bits).
   */
  input_size(): number {
    return this._m;
  }

  /**
   * Return the output size of this S-Box (in bits).
   */
  output_size(): number {
    return this._n;
  }

  /**
   * Apply the S-box to an input.
   *
   * If the input is a number, returns a number.
   * If the input is an array of bits, returns an array of bits.
   *
   * @param X - Either an integer or a tuple/array of bits
   * @returns The S-box output
   *
   * @example
   * ```typescript
   * const S = new SBox([7, 6, 0, 4, 2, 5, 1, 3]);
   * S.call(0);           // 7
   * S.call([0, 0, 0]);   // [1, 1, 1]
   * S.call([0, 0, 1]);   // [1, 1, 0]
   * ```
   */
  call(X: number | number[]): number | number[] {
    if (typeof X === 'number') {
      if (X < 0 || X >= this._S_list.length) {
        throw new ValueError(`input ${X} out of range [0, ${this._S_list.length - 1}]`);
      }
      return this._S_list[X]!;
    }

    // X is an array of bits
    if (X.length !== this._m) {
      throw new SageTypeError(`expected ${this._m} input bits, got ${X.length}`);
    }

    const index = this.from_bits(X);
    const output = this._S_list[index]!;
    return this.to_bits(output, this._n);
  }

  /**
   * Get the S-box output at index X (alias for call).
   */
  get(X: number): number {
    if (X < 0 || X >= this._S_list.length) {
      throw new ValueError(`index ${X} out of range [0, ${this._S_list.length - 1}]`);
    }
    return this._S_list[X]!;
  }

  /**
   * Return bitstring of length n for integer x.
   *
   * @param x - The integer to convert
   * @param n - The desired bit length (defaults to output size)
   *
   * @example
   * ```typescript
   * const S = new SBox([7, 6, 0, 4, 2, 5, 1, 3]);
   * S.to_bits(6);     // [1, 1, 0]
   * S.to_bits(6, 4);  // [0, 1, 1, 0]
   * ```
   */
  to_bits(x: number, n?: number): number[] {
    const bitLen = n ?? (this._m === this._n ? this._n : this._n);

    const bits: number[] = [];
    for (let i = 0; i < bitLen; i++) {
      bits.push((x >> i) & 1);
    }

    if (this._big_endian) {
      bits.reverse();
    }

    return bits;
  }

  /**
   * Return integer for bitstring x.
   *
   * @param x - A bitstring (array of 0s and 1s)
   * @param n - Expected bit length (optional)
   *
   * @example
   * ```typescript
   * const S = new SBox([7, 6, 0, 4, 2, 5, 1, 3]);
   * S.from_bits([1, 1, 0]);  // 6
   * ```
   */
  from_bits(x: number[], n?: number): number {
    const arr = this._big_endian ? [...x].reverse() : x;
    let result = 0;
    for (let i = 0; i < arr.length; i++) {
      if (arr[i]) {
        result |= 1 << i;
      }
    }
    return result;
  }

  /**
   * Return true if this S-Box is a permutation (bijective).
   *
   * @example
   * ```typescript
   * const S1 = new SBox([7, 6, 0, 4, 2, 5, 1, 3]);
   * S1.is_permutation();  // true
   *
   * const S2 = new SBox([3, 2, 0, 0, 2, 1, 1, 3]);
   * S2.is_permutation();  // false
   * ```
   */
  is_permutation(): boolean {
    if (this._m !== this._n) {
      return false;
    }

    const seen = new Set<number>();
    for (let i = 0; i < 1 << this._m; i++) {
      seen.add(this._S_list[i]!);
    }
    return seen.size === 1 << this._m;
  }

  /**
   * Return true if this S-Box is an involution (S(S(x)) = x for all x).
   *
   * @example
   * ```typescript
   * // AES S-box inverse is also an involution
   * const S = new SBox([...]);
   * S.is_involution();
   * ```
   */
  is_involution(): boolean {
    if (!this.is_permutation()) {
      return false;
    }

    for (let x = 0; x < 1 << this._m; x++) {
      if (this._S_list[this._S_list[x]!]! !== x) {
        return false;
      }
    }
    return true;
  }

  /**
   * Return the inverse of this S-Box.
   *
   * @throws {SageTypeError} If the S-Box is not a permutation
   *
   * @example
   * ```typescript
   * const S = new SBox([0, 1, 3, 6, 7, 4, 5, 2]);
   * const Sinv = S.inverse();
   * // Sinv(S(i)) === i for all i
   * ```
   */
  inverse(): SBox {
    if (!this.is_permutation()) {
      throw new SageTypeError('S-Box must be a permutation');
    }

    const inv: number[] = new Array(1 << this._m);
    for (let i = 0; i < 1 << this._m; i++) {
      inv[this._S_list[i]!] = i;
    }

    return new SBox(inv, { big_endian: this._big_endian });
  }

  /**
   * Return the derivative in direction u.
   *
   * The derivative of F in direction u is defined as x -> F(x) + F(x + u)
   * (where + is XOR in GF(2)).
   *
   * @param u - Either an integer or an array of bits
   *
   * @example
   * ```typescript
   * const S = new SBox([0, 1, 2, 3]);
   * S.derivative(1);  // SBox([1, 1, 1, 1])
   * ```
   */
  derivative(u: number | number[]): SBox {
    const v = typeof u === 'number' ? u : this.from_bits(u);

    const result: number[] = [];
    for (let x = 0; x < 1 << this._m; x++) {
      result.push(this._S_list[x]! ^ this._S_list[x ^ v]!);
    }

    return new SBox(result, { big_endian: this._big_endian });
  }

  /**
   * Return the Difference Distribution Table (DDT) for this S-box.
   *
   * The rows encode input differences Delta_I and columns encode output
   * differences Delta_O. Entry DDT[a][b] counts how many x satisfy:
   *   S(x) XOR S(x XOR a) = b
   *
   * The first row (a=0) always has DDT[0][0] = 2^m and all other entries 0.
   *
   * @returns A 2^m x 2^n matrix of integers
   *
   * @example
   * ```typescript
   * const S = new SBox([7, 6, 0, 4, 2, 5, 1, 3]);
   * const ddt = S.difference_distribution_table();
   * // ddt.get(0, 0) === 8 (first entry is always 2^m)
   * ```
   */
  difference_distribution_table(): Matrix<IntegerElement> {
    if (this._ddt_cache !== null) {
      return this._ddt_cache;
    }

    const nrows = 1 << this._m;
    const ncols = 1 << this._n;

    // Initialize the DDT with zeros
    const ddt: number[][] = [];
    for (let i = 0; i < nrows; i++) {
      ddt.push(new Array(ncols).fill(0));
    }

    // Compute DDT
    for (let x = 0; x < nrows; x++) {
      const sx = this._S_list[x]!;
      for (let a = 0; a < nrows; a++) {
        const diff_out = sx ^ this._S_list[x ^ a]!;
        ddt[a]![diff_out]!++;
      }
    }

    // Convert to Matrix
    const entries: IntegerElement[][] = [];
    for (let i = 0; i < nrows; i++) {
      entries.push(ddt[i]!.map((v) => new IntegerElement(BigInt(v))));
    }

    this._ddt_cache = new Matrix(IntegerRingForMatrix, nrows, ncols, entries);
    return this._ddt_cache;
  }

  /**
   * Return the differential uniformity of this S-Box.
   *
   * The differential uniformity is the maximum value in the DDT
   * excluding the trivial entry at (0, 0).
   *
   * A lower differential uniformity indicates better resistance to
   * differential cryptanalysis. The minimum possible value is 2.
   *
   * @example
   * ```typescript
   * const S = new SBox([7, 6, 0, 4, 2, 5, 1, 3]);
   * S.differential_uniformity();  // 2
   *
   * // AES S-box has differential uniformity 4
   * ```
   */
  differential_uniformity(): number {
    const ddt = this.difference_distribution_table();
    let maxVal = 0;

    const nrows = 1 << this._m;
    const ncols = 1 << this._n;

    for (let i = 0; i < nrows; i++) {
      for (let j = 0; j < ncols; j++) {
        // Skip the trivial entry (0, 0)
        if (i === 0 && j === 0) continue;
        const val = Number(ddt.get(i, j).value);
        if (val > maxVal) {
          maxVal = val;
        }
      }
    }

    return maxVal;
  }

  /**
   * Return true if this S-Box is Almost Perfect Nonlinear (APN).
   *
   * An m x m S-Box is APN if for every nonzero a and every b,
   * the equation S(x) XOR S(x XOR a) = b has at most 2 solutions.
   * Equivalently, the differential uniformity equals 2.
   *
   * @throws {SageTypeError} If input_size != output_size
   *
   * @example
   * ```typescript
   * const S = new SBox([0, 1, 3, 6, 7, 4, 5, 2]);
   * S.is_apn();  // true
   * S.differential_uniformity();  // 2
   * ```
   */
  is_apn(): boolean {
    if (this._m !== this._n) {
      throw new SageTypeError(
        'APN function is only defined for S-Boxes with equal input and output size'
      );
    }
    return this.differential_uniformity() === 2;
  }

  /**
   * Return the Linear Approximation Table (LAT) for this S-box.
   *
   * Entry LAT[a][b] is the bias of the linear approximation:
   *   LAT[a][b] = #{x : a.x = b.S(x)} - 2^{m-1}
   *
   * where a.x denotes the dot product in GF(2) (XOR of ANDed bits).
   *
   * The LAT is used in linear cryptanalysis. Entries close to 0 indicate
   * good resistance to linear cryptanalysis.
   *
   * @returns A 2^m x 2^n matrix of integers
   *
   * @example
   * ```typescript
   * const S = new SBox([7, 6, 0, 4, 2, 5, 1, 3]);
   * const lat = S.linear_approximation_table();
   * ```
   */
  linear_approximation_table(): Matrix<IntegerElement> {
    if (this._lat_cache !== null) {
      return this._lat_cache;
    }

    const nrows = 1 << this._m;
    const ncols = 1 << this._n;

    // Initialize LAT
    const lat: number[][] = [];
    for (let a = 0; a < nrows; a++) {
      lat.push(new Array(ncols).fill(0));
    }

    // Compute LAT using Walsh-Hadamard transform approach
    // LAT[a][b] = #{x : a.x = b.S(x)} - 2^{m-1}
    for (let a = 0; a < nrows; a++) {
      for (let b = 0; b < ncols; b++) {
        let count = 0;
        for (let x = 0; x < nrows; x++) {
          const ax = gf2DotProduct(a, x);
          const bSx = gf2DotProduct(b, this._S_list[x]!);
          if (ax === bSx) {
            count++;
          }
        }
        lat[a]![b] = count - (nrows >> 1);
      }
    }

    // Convert to Matrix
    const entries: IntegerElement[][] = [];
    for (let i = 0; i < nrows; i++) {
      entries.push(lat[i]!.map((v) => new IntegerElement(BigInt(v))));
    }

    this._lat_cache = new Matrix(IntegerRingForMatrix, nrows, ncols, entries);
    return this._lat_cache;
  }

  /**
   * Return the maximum absolute linear bias.
   *
   * This is the maximum absolute value in the LAT, excluding the
   * trivial entry at (0, 0).
   *
   * @example
   * ```typescript
   * const S = new SBox([7, 6, 0, 4, 2, 5, 1, 3]);
   * S.max_linear_bias();  // 2
   * ```
   */
  max_linear_bias(): number {
    const lat = this.linear_approximation_table();
    let maxAbs = 0;

    const nrows = 1 << this._m;
    const ncols = 1 << this._n;

    for (let i = 0; i < nrows; i++) {
      for (let j = 0; j < ncols; j++) {
        // Skip trivial entry
        if (i === 0 && j === 0) continue;
        const val = Math.abs(Number(lat.get(i, j).value));
        if (val > maxAbs) {
          maxAbs = val;
        }
      }
    }

    return maxAbs;
  }

  /**
   * Return the linearity of this S-Box.
   *
   * The linearity is 2 * max_linear_bias().
   *
   * @example
   * ```typescript
   * const S = new SBox([7, 6, 0, 4, 2, 5, 1, 3]);
   * S.linearity();  // 4
   * ```
   */
  linearity(): number {
    return 2 * this.max_linear_bias();
  }

  /**
   * Return the nonlinearity of this S-Box.
   *
   * The nonlinearity is defined as 2^{m-1} - max_linear_bias().
   *
   * Higher nonlinearity indicates better resistance to linear cryptanalysis.
   *
   * @example
   * ```typescript
   * const S = new SBox([7, 6, 0, 4, 2, 5, 1, 3]);
   * S.nonlinearity();  // 2
   * ```
   */
  nonlinearity(): number {
    return (1 << (this._m - 1)) - this.max_linear_bias();
  }

  /**
   * Return a Boolean function representing the component function b.S(x).
   *
   * The component function is the Boolean function f(x) = b.S(x),
   * where b.S(x) is the dot product of b and S(x) in GF(2).
   *
   * @param b - Either an integer or an array of bits specifying the output mask
   * @returns A BooleanFunction object representing the component function
   *
   * @example
   * ```typescript
   * const S = new SBox([7, 6, 0, 4, 2, 5, 1, 3]);
   * const f3 = S.component_function(3);  // b = 3 = 011 in binary
   * // f3 is the Boolean function x -> (S(x)[1] XOR S(x)[0])
   * console.log(f3.algebraicNormalForm());  // Shows ANF representation
   * ```
   */
  component_function(b: number | number[]): BooleanFunction {
    const mask = typeof b === 'number' ? b : this.from_bits(b, this._n);

    const truthTable: number[] = [];
    for (let x = 0; x < 1 << this._m; x++) {
      truthTable.push(gf2DotProduct(mask, this._S_list[x]!));
    }

    return new BooleanFunction(truthTable);
  }

  /**
   * Return the component function truth table as an array.
   *
   * This is a convenience method that returns the raw truth table instead
   * of a BooleanFunction object.
   *
   * @param b - Either an integer or an array of bits specifying the output mask
   * @returns An array of 0s and 1s representing the truth table
   */
  component_function_table(b: number | number[]): number[] {
    const mask = typeof b === 'number' ? b : this.from_bits(b, this._n);

    const result: number[] = [];
    for (let x = 0; x < 1 << this._m; x++) {
      result.push(gf2DotProduct(mask, this._S_list[x]!));
    }

    return result;
  }

  /**
   * Return the algebraic degree of the S-Box.
   *
   * The algebraic degree is the maximum degree of all component functions.
   *
   * For a Boolean function represented as a truth table, the algebraic degree
   * is computed from its Algebraic Normal Form (ANF).
   *
   * @example
   * ```typescript
   * const S = new SBox([12, 5, 6, 11, 9, 0, 10, 13, 3, 14, 15, 8, 4, 7, 1, 2]);
   * S.algebraic_degree();  // 3
   * ```
   */
  algebraic_degree(): number {
    let maxDeg = 0;

    // Check all non-zero component functions
    for (let b = 1; b < 1 << this._n; b++) {
      const f = this.component_function(b);
      const deg = f.algebraicDegree();
      if (deg > maxDeg) {
        maxDeg = deg;
      }
    }

    return maxDeg;
  }

  /**
   * Return the maximum algebraic degree of all component functions.
   *
   * Alias for `algebraic_degree()` for SageMath compatibility.
   *
   * @example
   * ```typescript
   * const S = new SBox([12, 5, 6, 11, 9, 0, 10, 13, 3, 14, 15, 8, 4, 7, 1, 2]);
   * S.max_degree();  // 3
   * ```
   */
  max_degree(): number {
    return this.algebraic_degree();
  }

  /**
   * Return the minimum algebraic degree of all component functions.
   *
   * Note: The zero function is considered to have degree -1 by convention,
   * but for S-box analysis we typically consider only non-constant component
   * functions. If all component functions are constant, returns 0.
   *
   * @example
   * ```typescript
   * const S = new SBox([12, 5, 6, 11, 9, 0, 10, 13, 3, 14, 15, 8, 4, 7, 1, 2]);
   * S.min_degree();  // 2
   * ```
   */
  min_degree(): number {
    let minDeg = this._m; // Maximum possible degree
    let foundNonConstant = false;

    // Check all non-zero component functions
    for (let b = 1; b < 1 << this._n; b++) {
      const f = this.component_function(b);
      const deg = f.algebraicDegree();

      // Skip constant functions (degree 0) or zero function (degree -1)
      // for finding the minimum of non-constant functions
      if (deg >= 1) {
        foundNonConstant = true;
        if (deg < minDeg) {
          minDeg = deg;
        }
      } else if (deg === 0) {
        // Constant non-zero function
        foundNonConstant = true;
        if (0 < minDeg) {
          minDeg = 0;
        }
      }
      // deg === -1 means zero function, skip it
    }

    // If all component functions are zero or constant, return 0
    return foundNonConstant ? minDeg : 0;
  }

  /**
   * Return the list of all fixed points of this S-Box.
   *
   * A fixed point is a value x such that S(x) = x.
   *
   * @example
   * ```typescript
   * const S = new SBox([0, 1, 3, 6, 7, 4, 5, 2]);
   * S.fixed_points();  // [0, 1]
   * ```
   */
  fixed_points(): number[] {
    const result: number[] = [];
    for (let x = 0; x < 1 << this._m; x++) {
      if (this._S_list[x] === x) {
        result.push(x);
      }
    }
    return result;
  }

  /**
   * Return the differential branch number of this S-Box.
   *
   * The differential branch number is:
   *   min_{v != w} { wt(v XOR w) + wt(S(v) XOR S(w)) }
   *
   * where wt(x) is the Hamming weight.
   */
  differential_branch_number(): number {
    let minBranch = (1 << this._m) + (1 << this._n);

    for (let a = 0; a < 1 << this._m; a++) {
      for (let b = 0; b < 1 << this._m; b++) {
        if (a !== b) {
          const inputDiff = a ^ b;
          const outputDiff = this._S_list[a]! ^ this._S_list[b]!;
          const branch = hammingWeight(inputDiff) + hammingWeight(outputDiff);
          if (branch < minBranch) {
            minBranch = branch;
          }
        }
      }
    }

    return minBranch;
  }

  /**
   * Return the linear branch number of this S-Box.
   *
   * The linear branch number is:
   *   min_{a != 0, b : LAT[a][b] != 0} { wt(a) + wt(b) }
   */
  linear_branch_number(): number {
    const lat = this.linear_approximation_table();
    let minBranch = (1 << this._m) + (1 << this._n);

    const nrows = 1 << this._m;
    const ncols = 1 << this._n;

    for (let a = 0; a < nrows; a++) {
      for (let b = 1; b < ncols; b++) {
        // Start b from 1 for non-trivial masks
        if (Number(lat.get(a, b).value) !== 0) {
          const branch = hammingWeight(a) + hammingWeight(b);
          if (branch < minBranch) {
            minBranch = branch;
          }
        }
      }
    }

    return minBranch;
  }

  /**
   * Check if this S-Box is balanced.
   *
   * An S-Box is balanced if all its component functions are balanced,
   * meaning each output value appears equally often.
   */
  is_balanced(): boolean {
    for (let b = 1; b < 1 << this._n; b++) {
      const f = this.component_function(b);
      if (!f.isBalanced()) {
        return false;
      }
    }
    return true;
  }

  /**
   * Check equality with another S-Box.
   */
  equals(other: SBox): boolean {
    if (this._m !== other._m || this._n !== other._n) {
      return false;
    }
    if (this._big_endian !== other._big_endian) {
      return false;
    }
    for (let i = 0; i < this._S_list.length; i++) {
      if (this._S_list[i] !== other._S_list[i]) {
        return false;
      }
    }
    return true;
  }

  /**
   * Iterate over all S-box outputs.
   */
  *[Symbol.iterator](): Iterator<number> {
    for (let i = 0; i < 1 << this._m; i++) {
      yield this._S_list[i]!;
    }
  }

  /**
   * Return a string representation of the S-Box.
   */
  toString(): string {
    return '(' + this._S_list.join(', ') + ')';
  }

  /**
   * Return the lookup table as an array.
   */
  toArray(): number[] {
    return [...this._S_list];
  }
}

// ============================================================================
// Standard S-Boxes
// ============================================================================

/**
 * The AES (Rijndael) S-Box.
 *
 * This is an 8x8 S-Box with differential uniformity 4 and
 * maximum linear bias 16 (linearity 32).
 *
 * The AES S-Box is constructed as the composition of:
 * 1. Multiplicative inverse in GF(2^8) (with 0 mapping to 0)
 * 2. An affine transformation over GF(2)
 */
export const AES_SBOX = new SBox([
  0x63, 0x7c, 0x77, 0x7b, 0xf2, 0x6b, 0x6f, 0xc5, 0x30, 0x01, 0x67, 0x2b, 0xfe, 0xd7, 0xab, 0x76,
  0xca, 0x82, 0xc9, 0x7d, 0xfa, 0x59, 0x47, 0xf0, 0xad, 0xd4, 0xa2, 0xaf, 0x9c, 0xa4, 0x72, 0xc0,
  0xb7, 0xfd, 0x93, 0x26, 0x36, 0x3f, 0xf7, 0xcc, 0x34, 0xa5, 0xe5, 0xf1, 0x71, 0xd8, 0x31, 0x15,
  0x04, 0xc7, 0x23, 0xc3, 0x18, 0x96, 0x05, 0x9a, 0x07, 0x12, 0x80, 0xe2, 0xeb, 0x27, 0xb2, 0x75,
  0x09, 0x83, 0x2c, 0x1a, 0x1b, 0x6e, 0x5a, 0xa0, 0x52, 0x3b, 0xd6, 0xb3, 0x29, 0xe3, 0x2f, 0x84,
  0x53, 0xd1, 0x00, 0xed, 0x20, 0xfc, 0xb1, 0x5b, 0x6a, 0xcb, 0xbe, 0x39, 0x4a, 0x4c, 0x58, 0xcf,
  0xd0, 0xef, 0xaa, 0xfb, 0x43, 0x4d, 0x33, 0x85, 0x45, 0xf9, 0x02, 0x7f, 0x50, 0x3c, 0x9f, 0xa8,
  0x51, 0xa3, 0x40, 0x8f, 0x92, 0x9d, 0x38, 0xf5, 0xbc, 0xb6, 0xda, 0x21, 0x10, 0xff, 0xf3, 0xd2,
  0xcd, 0x0c, 0x13, 0xec, 0x5f, 0x97, 0x44, 0x17, 0xc4, 0xa7, 0x7e, 0x3d, 0x64, 0x5d, 0x19, 0x73,
  0x60, 0x81, 0x4f, 0xdc, 0x22, 0x2a, 0x90, 0x88, 0x46, 0xee, 0xb8, 0x14, 0xde, 0x5e, 0x0b, 0xdb,
  0xe0, 0x32, 0x3a, 0x0a, 0x49, 0x06, 0x24, 0x5c, 0xc2, 0xd3, 0xac, 0x62, 0x91, 0x95, 0xe4, 0x79,
  0xe7, 0xc8, 0x37, 0x6d, 0x8d, 0xd5, 0x4e, 0xa9, 0x6c, 0x56, 0xf4, 0xea, 0x65, 0x7a, 0xae, 0x08,
  0xba, 0x78, 0x25, 0x2e, 0x1c, 0xa6, 0xb4, 0xc6, 0xe8, 0xdd, 0x74, 0x1f, 0x4b, 0xbd, 0x8b, 0x8a,
  0x70, 0x3e, 0xb5, 0x66, 0x48, 0x03, 0xf6, 0x0e, 0x61, 0x35, 0x57, 0xb9, 0x86, 0xc1, 0x1d, 0x9e,
  0xe1, 0xf8, 0x98, 0x11, 0x69, 0xd9, 0x8e, 0x94, 0x9b, 0x1e, 0x87, 0xe9, 0xce, 0x55, 0x28, 0xdf,
  0x8c, 0xa1, 0x89, 0x0d, 0xbf, 0xe6, 0x42, 0x68, 0x41, 0x99, 0x2d, 0x0f, 0xb0, 0x54, 0xbb, 0x16,
]);

/**
 * The PRESENT S-Box.
 *
 * This is a 4x4 S-Box used in the PRESENT lightweight block cipher.
 * It has differential uniformity 4 and maximum linear bias 4 (linearity 8).
 */
export const PRESENT_SBOX = new SBox([
  0xc, 0x5, 0x6, 0xb, 0x9, 0x0, 0xa, 0xd, 0x3, 0xe, 0xf, 0x8, 0x4, 0x7, 0x1, 0x2,
]);

/**
 * The DES (Data Encryption Standard) S-Box 1.
 *
 * This is a 6x4 S-Box (6 input bits, 4 output bits).
 */
export const DES_SBOX1 = new SBox([
  14, 4, 13, 1, 2, 15, 11, 8, 3, 10, 6, 12, 5, 9, 0, 7, 0, 15, 7, 4, 14, 2, 13, 1, 10, 6, 12, 11, 9,
  5, 3, 8, 4, 1, 14, 8, 13, 6, 2, 11, 15, 12, 9, 7, 3, 10, 5, 0, 15, 12, 8, 2, 4, 9, 1, 7, 5, 11, 3,
  14, 10, 0, 6, 13,
]);

/**
 * The GIFT S-Box.
 *
 * This is a 4x4 S-Box used in the GIFT lightweight block cipher.
 */
export const GIFT_SBOX = new SBox([
  0x1, 0xa, 0x4, 0xc, 0x6, 0xf, 0x3, 0x9, 0x2, 0xd, 0xb, 0x7, 0x5, 0x0, 0x8, 0xe,
]);

/**
 * The SKINNY S-Box.
 *
 * This is a 4x4 S-Box used in the SKINNY lightweight block cipher.
 */
export const SKINNY_SBOX = new SBox([
  0xc, 0x6, 0x9, 0x0, 0x1, 0xa, 0x2, 0xb, 0x3, 0x8, 0x5, 0xd, 0x4, 0xe, 0x7, 0xf,
]);

/**
 * Feistel construction for building larger S-boxes from smaller ones.
 *
 * Given S-boxes s1, s2, ..., sk (all with the same input/output size),
 * constructs an S-box with twice the input/output size using a k-round
 * Feistel network.
 *
 * @param sboxes - Array of S-boxes to use in the Feistel rounds
 * @returns A new S-box with doubled size
 *
 * @example
 * ```typescript
 * const s = new SBox([12, 5, 6, 11, 9, 0, 10, 13, 3, 14, 15, 8, 4, 7, 1, 2]);
 * const S = feistel_construction([s, s, s]);  // 8x8 S-box
 * ```
 */
export function feistel_construction(sboxes: SBox[]): SBox {
  if (sboxes.length === 0) {
    throw new SageTypeError('no input S-boxes provided');
  }

  const inputSize = sboxes[0]!.input_size();
  for (const sb of sboxes) {
    if (sb.input_size() !== inputSize || sb.output_size() !== inputSize) {
      throw new SageTypeError('all S-boxes must have the same input and output size');
    }
  }

  const mask = (1 << inputSize) - 1;
  const newSize = 2 * inputSize;
  const result: number[] = [];

  for (let x = 0; x < 1 << newSize; x++) {
    let left = (x >> inputSize) & mask;
    let right = x & mask;

    for (const sb of sboxes) {
      const newLeft = (sb.call(left) as number) ^ right;
      right = left;
      left = newLeft;
    }

    result.push((left << inputSize) | right);
  }

  return new SBox(result);
}

/**
 * Misty construction for building larger S-boxes from smaller ones.
 *
 * Similar to Feistel but with a different structure where the round
 * function is applied to the right half before XOR.
 *
 * @param sboxes - Array of S-boxes to use in the Misty rounds
 * @returns A new S-box with doubled size
 */
export function misty_construction(sboxes: SBox[]): SBox {
  if (sboxes.length === 0) {
    throw new SageTypeError('no input S-boxes provided');
  }

  const inputSize = sboxes[0]!.input_size();
  for (const sb of sboxes) {
    if (sb.input_size() !== inputSize || sb.output_size() !== inputSize) {
      throw new SageTypeError('all S-boxes must have the same input and output size');
    }
  }

  const mask = (1 << inputSize) - 1;
  const newSize = 2 * inputSize;
  const result: number[] = [];

  for (let x = 0; x < 1 << newSize; x++) {
    let left = (x >> inputSize) & mask;
    let right = x & mask;

    for (const sb of sboxes) {
      const newLeft = (sb.call(right) as number) ^ left;
      left = right;
      right = newLeft;
    }

    // Note: output is swapped compared to Feistel
    result.push((left << inputSize) | right);
  }

  return new SBox(result);
}
