/**
 * @module sage/crypto/boolean_function
 * @description Boolean functions for cryptographic analysis
 *
 * Boolean functions are used in cryptographic applications such as LFSR-based
 * ciphers (filter generators, combination generators). This module allows
 * studying spectral properties (Walsh-Hadamard transform) and algebraic
 * properties (algebraic immunity, algebraic degree).
 *
 * Port of: sage/crypto/boolean_function.pyx
 *
 * @example
 * ```typescript
 * import { BooleanFunction, randomBooleanFunction } from 'sagemath-ts/crypto/boolean_function';
 *
 * // Create from truth table
 * const f = new BooleanFunction([0, 1, 1, 0, 1, 0, 0, 1]);
 *
 * // Analyze cryptographic properties
 * console.log(f.nonlinearity());     // Distance to affine functions
 * console.log(f.isBalanced());       // True if equal 0s and 1s
 * console.log(f.algebraicDegree());  // Degree of ANF
 *
 * // Walsh-Hadamard transform for spectral analysis
 * const spectrum = f.walshHadamardTransform();
 * ```
 */

import { ValueError, NotImplementedError } from '../errors.js';
import { current_randstate } from '../misc/randstate.js';
import { GF2, GF2Element, GF2Field } from '../rings/finite_rings/gf2.js';

/**
 * Type for truth table input: array of 0/1 or boolean values.
 */
export type TruthTableInput = (0 | 1 | boolean | number)[];

/**
 * Represents a Boolean function f: GF(2)^n -> GF(2).
 *
 * A Boolean function on n variables can be represented by its truth table,
 * which has 2^n entries. It can also be represented by its Algebraic Normal
 * Form (ANF), which is a polynomial over GF(2).
 *
 * @example
 * ```typescript
 * // Create from truth table (XOR function on 2 variables)
 * const xor = new BooleanFunction([0, 1, 1, 0]);
 *
 * // Evaluate at specific inputs
 * console.log(xor.call(0));      // 0
 * console.log(xor.call(1));      // 1
 * console.log(xor.call(2));      // 1
 * console.log(xor.call(3));      // 0
 *
 * // Get Walsh transform
 * console.log(xor.walshHadamardTransform());  // [0, 0, 0, -4]
 * ```
 */
export class BooleanFunction {
  /**
   * The truth table stored as an array of booleans.
   * Index i represents f(i) where i is interpreted as an n-bit binary number.
   */
  private readonly _truthTable: boolean[];

  /**
   * Number of input variables n. The truth table has 2^n entries.
   */
  private readonly _nvariables: number;

  /**
   * Cached Walsh-Hadamard transform.
   */
  private _walshHadamardTransform: number[] | null = null;

  /**
   * Cached nonlinearity.
   */
  private _nonlinearity: number | null = null;

  /**
   * Cached correlation immunity order.
   */
  private _correlationImmunity: number | null = null;

  /**
   * Cached autocorrelation.
   */
  private _autocorrelation: number[] | null = null;

  /**
   * Cached algebraic degree.
   */
  private _algebraicDegree: number | null = null;

  /**
   * Construct a Boolean function.
   *
   * @param x - Can be:
   *   - An integer n: creates the zero function on n variables
   *   - An array of 0/1 or boolean values: the truth table (length must be 2^n)
   *   - A string in hexadecimal: represents the truth table
   *
   * @throws ValueError if truth table length is not a power of 2
   *
   * @example
   * ```typescript
   * // Zero function on 3 variables
   * const zero = new BooleanFunction(3);
   *
   * // From truth table
   * const f = new BooleanFunction([0, 1, 1, 0, 1, 0, 0, 1]);
   *
   * // From hex string (4 variables = 16 bits = 4 hex chars)
   * const g = new BooleanFunction('6996');
   * ```
   */
  constructor(x: number | TruthTableInput | string | BooleanFunction) {
    if (typeof x === 'number') {
      // Zero function with x variables
      if (x < 0 || !Number.isInteger(x)) {
        throw new ValueError('number of variables must be a non-negative integer');
      }
      this._nvariables = x;
      this._truthTable = new Array(1 << x).fill(false);
    } else if (typeof x === 'string') {
      // Hexadecimal truth table
      const len = x.length;
      if (len === 0 || (len & (len - 1)) !== 0) {
        throw new ValueError('the length of the truth table must be a power of 2');
      }

      // Each hex digit represents 4 bits
      const numBits = len * 4;
      this._nvariables = Math.log2(numBits);

      if (!Number.isInteger(this._nvariables)) {
        throw new ValueError('the length of the truth table must be a power of 2');
      }

      // Parse hex string to bits
      this._truthTable = [];
      for (let i = 0; i < len; i++) {
        const digit = parseInt(x[i]!, 16);
        if (isNaN(digit)) {
          throw new ValueError(`invalid hexadecimal character: ${x[i]}`);
        }
        // Extract 4 bits, LSB first within each nibble
        for (let j = 0; j < 4; j++) {
          this._truthTable.push(((digit >> j) & 1) === 1);
        }
      }
    } else if (x instanceof BooleanFunction) {
      // Copy constructor
      this._nvariables = x._nvariables;
      this._truthTable = [...x._truthTable];
    } else if (Array.isArray(x)) {
      // Truth table as array
      const len = x.length;
      if (len === 0 || (len & (len - 1)) !== 0) {
        throw new ValueError('the length of the truth table must be a power of 2');
      }

      this._nvariables = Math.log2(len);
      this._truthTable = x.map((v) => {
        if (typeof v === 'boolean') return v;
        return (Number(v) & 1) === 1;
      });
    } else {
      throw new ValueError('unable to init the Boolean function');
    }
  }

  /**
   * Clear cached values (called when the function is modified).
   */
  private _clearCache(): void {
    this._walshHadamardTransform = null;
    this._nonlinearity = null;
    this._correlationImmunity = null;
    this._autocorrelation = null;
    this._algebraicDegree = null;
  }

  /**
   * Return the number of input variables.
   *
   * @returns Number of variables n where input space is GF(2)^n
   *
   * @example
   * ```typescript
   * const f = new BooleanFunction([0, 1, 1, 0, 1, 0, 0, 1]);
   * console.log(f.nvariables());  // 3
   * ```
   */
  nvariables(): number {
    return this._nvariables;
  }

  /**
   * Return the size of the domain (2^n).
   *
   * @returns Number of possible input values
   */
  length(): number {
    return 1 << this._nvariables;
  }

  /**
   * Return the truth table as an array.
   *
   * @param format - Output format:
   *   - 'bool': array of booleans (default)
   *   - 'int': array of 0/1 integers
   *   - 'hex': hexadecimal string
   *
   * @returns Truth table in specified format
   *
   * @example
   * ```typescript
   * const f = new BooleanFunction([0, 1, 1, 0]);
   * console.log(f.truthTable());        // [false, true, true, false]
   * console.log(f.truthTable('int'));   // [0, 1, 1, 0]
   * console.log(f.truthTable('hex'));   // '6'
   * ```
   */
  truthTable(format: 'bool' | 'int' | 'hex' = 'bool'): boolean[] | number[] | string {
    if (format === 'bool') {
      return [...this._truthTable];
    }
    if (format === 'int') {
      return this._truthTable.map((b) => (b ? 1 : 0));
    }
    if (format === 'hex') {
      // Convert to hex string, LSB first within each nibble
      let result = '';
      const len = this._truthTable.length;
      for (let i = 0; i < len; i += 4) {
        let nibble = 0;
        for (let j = 0; j < 4 && i + j < len; j++) {
          if (this._truthTable[i + j]) {
            nibble |= 1 << j;
          }
        }
        result += nibble.toString(16);
      }
      return result;
    }
    throw new ValueError('unknown output format');
  }

  /**
   * Evaluate the function at a given input.
   *
   * @param x - Input value as integer or array of bits
   * @returns Function value (0 or 1)
   *
   * @throws ValueError if input is out of range
   *
   * @example
   * ```typescript
   * const f = new BooleanFunction([0, 1, 1, 0]);
   * console.log(f.call(0));        // 0
   * console.log(f.call(1));        // 1
   * console.log(f.call([0, 1]));   // 1 (same as call(2))
   * ```
   */
  call(x: number | number[]): 0 | 1 {
    let index: number;

    if (Array.isArray(x)) {
      if (x.length !== this._nvariables) {
        throw new ValueError(`expected ${this._nvariables} input bits, got ${x.length}`);
      }
      // Convert bit array to integer (LSB first)
      index = 0;
      for (let i = 0; i < x.length; i++) {
        if (x[i]) {
          index |= 1 << i;
        }
      }
    } else {
      index = x;
    }

    if (index < 0 || index >= this._truthTable.length) {
      throw new ValueError('index out of bound');
    }

    return this._truthTable[index] ? 1 : 0;
  }

  /**
   * Get the function value at index i.
   */
  get(i: number): boolean {
    return this._truthTable[i]!;
  }

  /**
   * Set the function value at index i.
   *
   * @param i - Input index
   * @param value - New value (0, 1, or boolean)
   */
  set(i: number, value: 0 | 1 | boolean): void {
    this._clearCache();
    this._truthTable[i] = typeof value === 'boolean' ? value : value === 1;
  }

  /**
   * Compute the Walsh-Hadamard transform of the function.
   *
   * The Walsh-Hadamard transform W_f is defined as:
   *   W_f(a) = sum_{x in GF(2)^n} (-1)^{f(x) + a.x}
   *
   * This uses the fast Walsh-Hadamard transform with O(n * 2^n) complexity.
   *
   * @returns Array of Walsh spectrum values W_f(0), W_f(1), ..., W_f(2^n - 1)
   *
   * @example
   * ```typescript
   * const f = new BooleanFunction([1, 0, 0, 1]);
   * const W = f.walshHadamardTransform();
   * console.log(W);  // [0, 0, 0, -4]
   *
   * // Verify Parseval's theorem: sum of W_f(a)^2 = 2^{2n}
   * const sumSquares = W.reduce((s, w) => s + w*w, 0);
   * console.log(sumSquares === 16);  // true (2^{2*2} = 16)
   * ```
   */
  walshHadamardTransform(): number[] {
    if (this._walshHadamardTransform !== null) {
      return [...this._walshHadamardTransform];
    }

    const n = this._truthTable.length;

    // Initialize: transform f(x) to (-1)^{f(x)}
    const f = new Array<number>(n);
    for (let i = 0; i < n; i++) {
      f[i] = this._truthTable[i] ? -1 : 1;
    }

    // Fast Walsh-Hadamard transform
    // At each stage ldm, we process blocks of size 2^ldm
    const ldn = this._nvariables;
    for (let ldm = 1; ldm <= ldn; ldm++) {
      const m = 1 << ldm;
      const mh = m >> 1;

      for (let r = 0; r < n; r += m) {
        for (let j = 0; j < mh; j++) {
          const t1 = r + j;
          const t2 = r + j + mh;
          const u = f[t1]!;
          const v = f[t2]!;
          f[t1] = u + v;
          f[t2] = u - v;
        }
      }
    }

    this._walshHadamardTransform = f;
    return [...f];
  }

  /**
   * Return the absolute Walsh spectrum as a frequency map.
   *
   * @returns Map from absolute Walsh values to their frequency
   *
   * @example
   * ```typescript
   * const f = new BooleanFunction('7969817CC5893BA6AC326E47619F5AD0');
   * const spectrum = f.absoluteWalshSpectrum();
   * console.log(spectrum);  // Map { 0 => 64, 16 => 64 }
   * ```
   */
  absoluteWalshSpectrum(): Map<number, number> {
    const W = this.walshHadamardTransform();
    const d = new Map<number, number>();

    for (const w of W) {
      const a = Math.abs(w);
      d.set(a, (d.get(a) ?? 0) + 1);
    }

    return d;
  }

  /**
   * Check if the function is balanced.
   *
   * A Boolean function is balanced if it takes the value 1 exactly half
   * the time (equivalently, if W_f(0) = 0).
   *
   * @returns True if the function is balanced
   *
   * @example
   * ```typescript
   * const f = new BooleanFunction([0, 1, 1, 0]);
   * console.log(f.isBalanced());  // true
   *
   * const g = new BooleanFunction([0, 0, 1, 1, 1, 1, 1, 1]);
   * console.log(g.isBalanced());  // false (6 ones, 2 zeros)
   * ```
   */
  isBalanced(): boolean {
    // A function is balanced iff W_f(0) = 0
    return this.walshHadamardTransform()[0] === 0;
  }

  /**
   * Compute the nonlinearity of the function.
   *
   * The nonlinearity is the distance to the nearest affine function,
   * computed as: NL(f) = 2^{n-1} - (1/2) * max_a |W_f(a)|
   *
   * For optimal security against linear attacks, higher nonlinearity is better.
   * The maximum possible nonlinearity for n variables is 2^{n-1} - 2^{n/2-1}
   * (achieved by bent functions when n is even).
   *
   * @returns The nonlinearity value
   *
   * @example
   * ```typescript
   * // Bent function on 4 variables has maximum nonlinearity 6
   * const bent = new BooleanFunction('0113077C165E76A8');
   * console.log(bent.nonlinearity());  // 28 (max for 8 variables)
   * ```
   */
  nonlinearity(): number {
    if (this._nonlinearity !== null) {
      return this._nonlinearity;
    }

    const W = this.walshHadamardTransform();
    let maxAbs = 0;
    for (const w of W) {
      const abs = Math.abs(w);
      if (abs > maxAbs) maxAbs = abs;
    }

    this._nonlinearity = ((1 << this._nvariables) - maxAbs) >> 1;
    return this._nonlinearity;
  }

  /**
   * Check if the function is bent.
   *
   * A bent function has the maximum possible nonlinearity for its number of
   * variables. Bent functions exist only for even numbers of variables.
   *
   * @returns True if the function is bent
   *
   * @example
   * ```typescript
   * const f = new BooleanFunction('0113077C165E76A8');
   * console.log(f.isBent());  // true
   * ```
   */
  isBent(): boolean {
    // Bent functions only exist for even n
    if (this._nvariables & 1) {
      return false;
    }
    // Bent function has nonlinearity 2^{n-1} - 2^{n/2-1}
    const maxNL = (1 << (this._nvariables - 1)) - (1 << (this._nvariables / 2 - 1));
    return this.nonlinearity() === maxNL;
  }

  /**
   * Compute the correlation immunity order of the function.
   *
   * A Boolean function is correlation immune of order m if its output is
   * statistically independent of any combination of m input bits.
   * Equivalently, W_f(a) = 0 for all a with Hamming weight 1 <= wt(a) <= m.
   *
   * @returns The maximum order m of correlation immunity
   *
   * @example
   * ```typescript
   * const f = new BooleanFunction('7969817CC5893BA6AC326E47619F5AD0');
   * console.log(f.correlationImmunity());  // 2
   * ```
   */
  correlationImmunity(): number {
    if (this._correlationImmunity !== null) {
      return this._correlationImmunity;
    }

    const W = this.walshHadamardTransform();
    let minWeight = this._nvariables;

    for (let i = 1; i < W.length; i++) {
      if (W[i] !== 0) {
        const weight = hammingWeight(i);
        if (weight < minWeight) {
          minWeight = weight;
        }
      }
    }

    this._correlationImmunity = minWeight - 1;
    return this._correlationImmunity;
  }

  /**
   * Return the resiliency order of the function.
   *
   * A function is resilient of order m if it is balanced and correlation
   * immune of order m. If the function is not balanced, returns -1.
   *
   * @returns Resiliency order, or -1 if not balanced
   *
   * @example
   * ```typescript
   * const f = new BooleanFunction('077CE5A2F8831A5DF8831A5D077CE5A26996699669699696669999665AA5A55A');
   * console.log(f.resiliencyOrder());  // 3
   * ```
   */
  resiliencyOrder(): number {
    if (!this.isBalanced()) {
      return -1;
    }
    return this.correlationImmunity();
  }

  /**
   * Compute the autocorrelation of the function.
   *
   * The autocorrelation is defined as:
   *   Delta_f(a) = sum_{x in GF(2)^n} (-1)^{f(x) + f(x+a)}
   *
   * This measures the correlation between f(x) and f(x+a).
   *
   * @returns Array of autocorrelation values
   *
   * @example
   * ```typescript
   * const f = new BooleanFunction('03');
   * const ac = f.autocorrelation();
   * console.log(ac);  // [8, 8, 0, 0, 0, 0, 0, 0]
   * ```
   */
  autocorrelation(): number[] {
    if (this._autocorrelation !== null) {
      return [...this._autocorrelation];
    }

    const n = this._truthTable.length;
    const W = this.walshHadamardTransform();

    // Autocorrelation is the inverse Walsh transform of W^2
    const temp = new Array<number>(n);
    for (let i = 0; i < n; i++) {
      temp[i] = W[i]! * W[i]!;
    }

    // Apply Walsh-Hadamard transform again
    const ldn = this._nvariables;
    for (let ldm = 1; ldm <= ldn; ldm++) {
      const m = 1 << ldm;
      const mh = m >> 1;

      for (let r = 0; r < n; r += m) {
        for (let j = 0; j < mh; j++) {
          const t1 = r + j;
          const t2 = r + j + mh;
          const u = temp[t1]!;
          const v = temp[t2]!;
          temp[t1] = u + v;
          temp[t2] = u - v;
        }
      }
    }

    // Normalize by 2^n
    for (let i = 0; i < n; i++) {
      temp[i] = temp[i]! >> this._nvariables;
    }

    this._autocorrelation = temp;
    return [...temp];
  }

  /**
   * Return the absolute autocorrelation spectrum as a frequency map.
   *
   * @returns Map from absolute autocorrelation values to their frequency
   */
  absoluteAutocorrelation(): Map<number, number> {
    const ac = this.autocorrelation();
    const d = new Map<number, number>();

    for (const a of ac) {
      const absVal = Math.abs(a);
      d.set(absVal, (d.get(absVal) ?? 0) + 1);
    }

    return d;
  }

  /**
   * Return the absolute indicator (maximum absolute autocorrelation for a != 0).
   *
   * @returns Maximum |Delta_f(a)| for a != 0
   */
  absoluteIndicator(): number {
    const ac = this.autocorrelation();
    let maxAbs = 0;
    // Skip a=0 (ac[0] is always 2^n)
    for (let i = 1; i < ac.length; i++) {
      const abs = Math.abs(ac[i]!);
      if (abs > maxAbs) maxAbs = abs;
    }
    return maxAbs;
  }

  /**
   * Return the sum of square indicator.
   *
   * @returns Sum of Delta_f(a)^2 for all a
   */
  sumOfSquareIndicator(): number {
    const ac = this.autocorrelation();
    let sum = 0;
    for (const a of ac) {
      sum += a * a;
    }
    return sum;
  }

  /**
   * Check if the function is plateaued.
   *
   * A function is plateaued if its Walsh transform takes at most three
   * values: 0 and +/- lambda for some positive integer lambda.
   *
   * @returns True if the function is plateaued
   */
  isPlateaued(): boolean {
    const spectrum = this.absoluteWalshSpectrum();
    if (spectrum.size === 1) return true;
    if (spectrum.size === 2 && spectrum.has(0)) return true;
    return false;
  }

  /**
   * Compute the derivative of f in direction u.
   *
   * The derivative D_u f is defined as:
   *   D_u f(x) = f(x) + f(x + u)  (addition in GF(2)^n and GF(2))
   *
   * @param u - Direction vector as integer or bit array
   * @returns The derivative function
   *
   * @example
   * ```typescript
   * const f = new BooleanFunction([0, 1, 0, 1, 0, 1, 0, 1]);
   * const df = f.derivative(1);
   * console.log(df.truthTable('int'));  // [1, 1, 1, 1, 1, 1, 1, 1]
   * ```
   */
  derivative(u: number | number[]): BooleanFunction {
    let direction: number;

    if (Array.isArray(u)) {
      if (u.length !== this._nvariables) {
        throw new ValueError(`expected ${this._nvariables} bits, got ${u.length}`);
      }
      direction = 0;
      for (let i = 0; i < u.length; i++) {
        if (u[i]) {
          direction |= 1 << i;
        }
      }
    } else {
      direction = u;
    }

    if (direction < 0 || direction >= this._truthTable.length) {
      throw new ValueError('index out of bound');
    }

    const n = this._truthTable.length;
    const result: number[] = [];
    for (let x = 0; x < n; x++) {
      const fx = this._truthTable[x] ? 1 : 0;
      const fxu = this._truthTable[x ^ direction] ? 1 : 0;
      result.push(fx ^ fxu);
    }

    return new BooleanFunction(result);
  }

  /**
   * Check if a value is a linear structure of this function.
   *
   * A vector a is a linear structure if D_a f(x) is constant.
   *
   * @param val - Value to check (integer or bit array)
   * @returns True if val is a linear structure
   */
  isLinearStructure(val: number | number[]): boolean {
    let index: number;

    if (Array.isArray(val)) {
      if (val.length !== this._nvariables) {
        throw new ValueError(`expected ${this._nvariables} bits, got ${val.length}`);
      }
      index = 0;
      for (let i = 0; i < val.length; i++) {
        if (val[i]) {
          index |= 1 << i;
        }
      }
    } else {
      index = val;
    }

    if (index < 0 || index >= this._truthTable.length) {
      throw new ValueError('index out of range');
    }

    const ac = this.autocorrelation();
    return Math.abs(ac[index]!) === 1 << this._nvariables;
  }

  /**
   * Check if this function has any linear structure.
   *
   * @returns True if there exists a nonzero linear structure
   */
  hasLinearStructure(): boolean {
    const ac = this.autocorrelation();
    const threshold = 1 << this._nvariables;

    for (let i = 1; i < ac.length; i++) {
      if (Math.abs(ac[i]!) === threshold) {
        return true;
      }
    }
    return false;
  }

  /**
   * Compute the Algebraic Normal Form (ANF) coefficients.
   *
   * The ANF is a polynomial representation over GF(2):
   *   f(x1,...,xn) = XOR of a_I * prod_{i in I} x_i
   *
   * This uses the Reed-Muller (binary Mobius) transform.
   *
   * @returns Array of ANF coefficients (coefficient at index i corresponds to monomial
   *          where x_j appears iff bit j of i is set)
   *
   * @example
   * ```typescript
   * const f = new BooleanFunction([0, 1, 1, 0]);
   * const anf = f.algebraicNormalFormCoefficients();
   * console.log(anf);  // [0, 1, 1, 0] means f = x0 + x1 (XOR)
   * ```
   */
  algebraicNormalFormCoefficients(): (0 | 1)[] {
    const n = this._truthTable.length;
    const anf: (0 | 1)[] = this._truthTable.map((b) => (b ? 1 : 0));

    // Reed-Muller transform (binary Mobius transform)
    for (let i = 0; i < this._nvariables; i++) {
      const step = 1 << i;
      for (let j = 0; j < n; j += 2 * step) {
        for (let k = 0; k < step; k++) {
          anf[j + k + step] = (anf[j + k + step]! ^ anf[j + k]!) as 0 | 1;
        }
      }
    }

    return anf;
  }

  /**
   * Return the algebraic normal form as a human-readable string.
   *
   * @param varNames - Variable names to use (default: x0, x1, ...)
   * @returns String representation of the ANF polynomial
   *
   * @example
   * ```typescript
   * const f = new BooleanFunction([0, 1, 1, 0, 1, 0, 1, 1]);
   * console.log(f.algebraicNormalForm());
   * // Output: x0 + x1*x2 + x0*x1*x2 + x1 + x2
   * ```
   */
  algebraicNormalForm(varNames?: string[]): string {
    const anf = this.algebraicNormalFormCoefficients();
    const names = varNames ?? Array.from({ length: this._nvariables }, (_, i) => `x${i}`);

    if (names.length !== this._nvariables) {
      throw new ValueError(`expected ${this._nvariables} variable names`);
    }

    const terms: string[] = [];

    for (let i = 0; i < anf.length; i++) {
      if (anf[i]) {
        if (i === 0) {
          terms.push('1');
        } else {
          const vars: string[] = [];
          for (let j = 0; j < this._nvariables; j++) {
            if ((i >> j) & 1) {
              vars.push(names[j]!);
            }
          }
          terms.push(vars.join('*'));
        }
      }
    }

    if (terms.length === 0) {
      return '0';
    }

    return terms.join(' + ');
  }

  /**
   * Compute the algebraic degree of the function.
   *
   * The algebraic degree is the maximum weight of a monomial in the ANF.
   * By convention, the zero function has degree -1.
   *
   * @returns Algebraic degree, or -1 for the zero function
   *
   * @example
   * ```typescript
   * const f = new BooleanFunction([0, 1, 1, 0, 1, 0, 1, 1]);
   * console.log(f.algebraicDegree());  // 3 (due to x0*x1*x2 term)
   * ```
   */
  algebraicDegree(): number {
    if (this._algebraicDegree !== null) {
      return this._algebraicDegree;
    }

    const anf = this.algebraicNormalFormCoefficients();
    let maxDegree = -1;

    for (let i = 0; i < anf.length; i++) {
      if (anf[i]) {
        const deg = hammingWeight(i);
        if (deg > maxDegree) {
          maxDegree = deg;
        }
      }
    }

    this._algebraicDegree = maxDegree;
    return maxDegree;
  }

  /**
   * Return the complement of this Boolean function.
   *
   * @returns New function g where g(x) = 1 - f(x)
   *
   * @example
   * ```typescript
   * const f = new BooleanFunction([0, 1, 1, 0]);
   * const g = f.complement();
   * console.log(g.truthTable('int'));  // [1, 0, 0, 1]
   * ```
   */
  complement(): BooleanFunction {
    const result = new BooleanFunction(this._nvariables);
    for (let i = 0; i < this._truthTable.length; i++) {
      result._truthTable[i] = !this._truthTable[i];
    }
    return result;
  }

  /**
   * Return the XOR (addition in GF(2)) of this function with another.
   *
   * @param other - Another Boolean function with the same number of variables
   * @returns New function f + g
   */
  add(other: BooleanFunction): BooleanFunction {
    if (this._nvariables !== other._nvariables) {
      throw new ValueError('Boolean functions must have the same number of variables');
    }

    const result = new BooleanFunction(this._nvariables);
    for (let i = 0; i < this._truthTable.length; i++) {
      result._truthTable[i] = this._truthTable[i] !== other._truthTable[i];
    }
    return result;
  }

  /**
   * Return the AND (multiplication in GF(2)) of this function with another.
   *
   * @param other - Another Boolean function with the same number of variables
   * @returns New function f * g
   */
  mul(other: BooleanFunction): BooleanFunction {
    if (this._nvariables !== other._nvariables) {
      throw new ValueError('Boolean functions must have the same number of variables');
    }

    const result = new BooleanFunction(this._nvariables);
    for (let i = 0; i < this._truthTable.length; i++) {
      result._truthTable[i] = this._truthTable[i] && other._truthTable[i];
    }
    return result;
  }

  /**
   * Return the concatenation of this function with another.
   *
   * Creates a new function on n+1 variables where the first half of the
   * truth table comes from this function and the second half from other.
   *
   * @param other - Another Boolean function with the same number of variables
   * @returns Concatenated function on n+1 variables
   */
  concatenate(other: BooleanFunction): BooleanFunction {
    if (this._nvariables !== other._nvariables) {
      throw new ValueError('Boolean functions must have the same number of variables');
    }

    const newTT = [...this._truthTable, ...other._truthTable];
    return new BooleanFunction(newTT);
  }

  /**
   * Check equality with another Boolean function.
   *
   * @param other - Another Boolean function
   * @returns True if the functions are identical
   */
  equals(other: BooleanFunction): boolean {
    if (this._nvariables !== other._nvariables) {
      return false;
    }
    for (let i = 0; i < this._truthTable.length; i++) {
      if (this._truthTable[i] !== other._truthTable[i]) {
        return false;
      }
    }
    return true;
  }

  /**
   * Check if this function is symmetric.
   *
   * A symmetric function's output depends only on the Hamming weight of
   * the input, not on the specific bit positions.
   *
   * @returns True if the function is symmetric
   */
  isSymmetric(): boolean {
    // For symmetric functions, f(x) depends only on wt(x)
    // Check that all inputs with the same weight give the same output
    const n = this._truthTable.length;
    const valuesByWeight = new Map<number, boolean | null>();

    for (let i = 0; i < n; i++) {
      const wt = hammingWeight(i);
      const val = this._truthTable[i]!;

      const existing = valuesByWeight.get(wt);
      if (existing === null) {
        // Already found inconsistency
        continue;
      } else if (existing === undefined) {
        valuesByWeight.set(wt, val);
      } else if (existing !== val) {
        return false;
      }
    }

    return true;
  }

  /**
   * Compute the algebraic immunity of the function.
   *
   * The algebraic immunity is the minimum degree d such that there exists
   * a nonzero function g of degree d with f*g = 0 or (1-f)*g = 0.
   *
   * This is a simplified implementation that may be slow for large n.
   *
   * @returns Algebraic immunity value
   */
  algebraicImmunity(): number {
    // Find minimum degree d such that f or (1-f) has an annihilator of degree d
    const fc = this.complement();

    for (let d = 0; d <= this._nvariables; d++) {
      if (this._hasAnnihilatorOfDegree(d) || fc._hasAnnihilatorOfDegree(d)) {
        return d;
      }
    }

    // Should never reach here
    return this._nvariables;
  }

  /**
   * Check if there exists a nonzero annihilator of given degree.
   * Internal helper for algebraicImmunity.
   */
  private _hasAnnihilatorOfDegree(d: number): boolean {
    // An annihilator g of degree <= d satisfies f*g = 0
    // This means g(x) = 0 whenever f(x) = 1

    // Get the support of f (positions where f = 1)
    const support: number[] = [];
    for (let i = 0; i < this._truthTable.length; i++) {
      if (this._truthTable[i]) {
        support.push(i);
      }
    }

    if (support.length === 0) {
      // Zero function has annihilator of degree 0 (the constant 1)
      return d >= 0;
    }

    // Generate all monomials of degree <= d
    const monomials: number[] = [];
    for (let i = 0; i < this._truthTable.length; i++) {
      if (hammingWeight(i) <= d) {
        monomials.push(i);
      }
    }

    if (monomials.length === 0) {
      return false;
    }

    // Build the constraint matrix: rows = support points, cols = monomials
    // We need to find a nonzero solution to M*x = 0 (mod 2)
    // where M[i][j] = monomial j evaluated at support[i]

    // Use Gaussian elimination over GF(2)
    const rows = support.length;
    const cols = monomials.length;

    if (cols <= rows) {
      // Not enough monomials to have a nonzero kernel
      return false;
    }

    // Build matrix
    const M: number[][] = [];
    for (let i = 0; i < rows; i++) {
      const row: number[] = [];
      const x = support[i]!;
      for (let j = 0; j < cols; j++) {
        // Evaluate monomial j at x
        // Monomial j is the product of x_k where bit k of monomials[j] is set
        const mono = monomials[j]!;
        // Value is 1 iff all variables in the monomial are 1 in x
        // i.e., (x & mono) === mono
        row.push((x & mono) === mono ? 1 : 0);
      }
      M.push(row);
    }

    // Gaussian elimination to find rank
    let rank = 0;
    for (let col = 0; col < cols && rank < rows; col++) {
      // Find pivot
      let pivot = -1;
      for (let row = rank; row < rows; row++) {
        if (M[row]![col]) {
          pivot = row;
          break;
        }
      }

      if (pivot === -1) continue;

      // Swap rows
      [M[rank], M[pivot]] = [M[pivot]!, M[rank]!];

      // Eliminate
      for (let row = 0; row < rows; row++) {
        if (row !== rank && M[row]![col]) {
          for (let c = 0; c < cols; c++) {
            M[row]![c] = M[row]![c]! ^ M[rank]![c]!;
          }
        }
      }

      rank++;
    }

    // Kernel dimension is cols - rank
    return cols > rank;
  }

  /**
   * String representation of the Boolean function.
   */
  toString(): string {
    const varText = this._nvariables === 1 ? 'variable' : 'variables';
    return `Boolean function with ${this._nvariables} ${varText}`;
  }

  /**
   * Iterate over truth table values.
   */
  *[Symbol.iterator](): Iterator<boolean> {
    for (const val of this._truthTable) {
      yield val;
    }
  }
}

/**
 * Compute the Hamming weight (number of 1 bits) of a non-negative integer.
 *
 * @param n - Non-negative integer
 * @returns Number of 1 bits in binary representation
 */
export function hammingWeight(n: number): number {
  let count = 0;
  while (n) {
    count += n & 1;
    n >>>= 1;
  }
  return count;
}

/**
 * Generate a random Boolean function on n variables.
 *
 * @param n - Number of variables
 * @returns Random Boolean function
 *
 * @example
 * ```typescript
 * const f = randomBooleanFunction(8);
 * console.log(f.nvariables());  // 8
 * console.log(f.nonlinearity());  // Typically around 100-110 for random 8-var functions
 * ```
 */
export function randomBooleanFunction(n: number): BooleanFunction {
  const len = 1 << n;
  const tt: number[] = [];
  for (let i = 0; i < len; i++) {
    tt.push(current_randstate().random() < 0.5 ? 0 : 1);
  }
  return new BooleanFunction(tt);
}

/**
 * Construct a Boolean function from its Algebraic Normal Form coefficients.
 *
 * @param coefficients - ANF coefficients (coefficient at index i corresponds to
 *                       monomial where x_j appears iff bit j of i is set)
 * @returns Boolean function with the given ANF
 *
 * @example
 * ```typescript
 * // Create f = x0 + x1 (XOR function)
 * const f = fromANF([0, 1, 1, 0]);
 * console.log(f.truthTable('int'));  // [0, 1, 1, 0]
 * ```
 */
export function fromANF(coefficients: (0 | 1 | number)[]): BooleanFunction {
  const len = coefficients.length;
  if (len === 0 || (len & (len - 1)) !== 0) {
    throw new ValueError('coefficient array length must be a power of 2');
  }

  const nvariables = Math.log2(len);

  // Convert ANF to truth table using inverse Reed-Muller transform
  // (same transform, it's self-inverse)
  const tt = coefficients.map((c) => (Number(c) & 1));

  for (let i = 0; i < nvariables; i++) {
    const step = 1 << i;
    for (let j = 0; j < len; j += 2 * step) {
      for (let k = 0; k < step; k++) {
        tt[j + k + step] = tt[j + k + step]! ^ tt[j + k]!;
      }
    }
  }

  return new BooleanFunction(tt);
}

/**
 * Create a bent function on n variables (n must be even).
 *
 * This creates the "simplest" bent function: the inner product x1*y1 + ... + xk*yk
 * where n = 2k.
 *
 * @param n - Number of variables (must be even and >= 2)
 * @returns A bent function
 *
 * @example
 * ```typescript
 * const f = createBentFunction(4);
 * console.log(f.isBent());  // true
 * console.log(f.nonlinearity());  // 6 (maximum for 4 variables)
 * ```
 */
export function createBentFunction(n: number): BooleanFunction {
  if (n < 2 || n % 2 !== 0) {
    throw new ValueError('bent functions require an even number of variables >= 2');
  }

  const k = n / 2;
  const len = 1 << n;
  const tt: number[] = [];

  for (let i = 0; i < len; i++) {
    // Split i into two k-bit parts
    const x = i & ((1 << k) - 1);
    const y = i >> k;
    // Compute inner product x.y
    let innerProduct = 0;
    for (let j = 0; j < k; j++) {
      innerProduct ^= ((x >> j) & 1) & ((y >> j) & 1);
    }
    tt.push(innerProduct);
  }

  return new BooleanFunction(tt);
}

/**
 * Create an affine function a.x + b.
 *
 * @param a - Linear part as integer (bit i indicates coefficient of x_i)
 * @param b - Constant term (0 or 1)
 * @param n - Number of variables
 * @returns Affine Boolean function
 *
 * @example
 * ```typescript
 * // Create f(x0,x1,x2) = x0 + x2 + 1
 * const f = createAffineFunction(0b101, 1, 3);
 * console.log(f.truthTable('int'));  // [1, 0, 1, 0, 0, 1, 0, 1]
 * ```
 */
export function createAffineFunction(a: number, b: 0 | 1, n: number): BooleanFunction {
  const len = 1 << n;
  const tt: number[] = [];

  for (let x = 0; x < len; x++) {
    // Compute a.x (inner product)
    const ax = hammingWeight(a & x) & 1;
    tt.push(ax ^ b);
  }

  return new BooleanFunction(tt);
}

/**
 * Verify Parseval's theorem for a Walsh-Hadamard transform.
 *
 * Parseval's theorem states that sum of W_f(a)^2 = 2^{2n} for any Boolean
 * function f on n variables.
 *
 * @param W - Walsh-Hadamard transform array
 * @returns True if Parseval's theorem holds
 *
 * @example
 * ```typescript
 * const f = new BooleanFunction([0, 1, 1, 0, 1, 0, 0, 1]);
 * const W = f.walshHadamardTransform();
 * console.log(verifyParseval(W));  // true
 * ```
 */
export function verifyParseval(W: number[]): boolean {
  const n = Math.log2(W.length);
  if (!Number.isInteger(n)) return false;

  let sumSquares = 0;
  for (const w of W) {
    sumSquares += w * w;
  }

  const expected = 1 << (2 * n); // 2^{2n}
  return sumSquares === expected;
}
