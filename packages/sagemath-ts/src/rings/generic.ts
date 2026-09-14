/**
 * Integer-domain product/remainder trees and streaming product-rule evaluation.
 * @see Reference: sage/rings/generic.py
 * @see Deviation: Product trees and factor-base nontermination
 */
import { CRT_basis } from '../arith/misc.js';
import { AssertionError, ValueError, ZeroDivisionError } from '../errors.js';
import { type IntegerLike, toBigInt } from '../types/coercion.js';

function remainder(value: bigint, modulus: bigint): bigint {
  if (modulus === 0n) throw new ZeroDivisionError('Integer modulo by zero');
  const r = value % modulus;
  return r !== 0n && r < 0n !== modulus < 0n ? r + modulus : r;
}

export class ProductTree implements Iterable<bigint> {
  layers: Array<readonly bigint[]>;
  _crt_bases: bigint[][][] | null = null;

  constructor(leaves: Iterable<IntegerLike>) {
    let layer: readonly bigint[] = Object.freeze(Array.from(leaves, toBigInt));
    this.layers = [layer];
    while (layer.length > 1) {
      const next: bigint[] = [];
      for (let i = 0; i < layer.length; i += 2) {
        next.push(i + 1 < layer.length ? layer[i]! * layer[i + 1]! : layer[i]!);
      }
      layer = Object.freeze(next);
      this.layers.push(layer);
    }
  }

  __len__(): number {
    return this.layers[0]!.length;
  }
  __iter__(): IterableIterator<bigint> {
    return this.layers[0]!.values();
  }
  [Symbol.iterator](): IterableIterator<bigint> {
    return this.__iter__();
  }

  root(): bigint {
    const layer = this.layers[this.layers.length - 1]!;
    if (layer.length !== 1) throw new AssertionError('');
    return layer[0]!;
  }
  leaves(): readonly bigint[] {
    return this.layers[0]!;
  }

  remainders(x: IntegerLike): bigint[] {
    let values = [toBigInt(x)];
    for (let j = this.layers.length - 1; j >= 0; j--) {
      values = this.layers[j]!.map((m, i) => remainder(values[Math.floor(i / 2)]!, m));
    }
    return values;
  }

  interpolation(xs: readonly IntegerLike[]): bigint {
    // Initialization precedes the length check. A failed initialization leaves
    // the same partial cache as the original, observable on subsequent calls.
    if (this._crt_bases === null) {
      this._crt_bases = [];
      for (const layer of this.layers.slice(0, -1)) {
        const basis: bigint[][] = [];
        for (let i = 0; i < layer.length; i += 2) {
          basis.push(CRT_basis(layer.slice(i, i + 2)) as bigint[]);
        }
        this._crt_bases.push(basis);
      }
    }
    if (xs.length !== this.layers[0]!.length) {
      throw new ValueError('number of given elements must equal the number of leaves');
    }
    let values = xs.map(toBigInt);
    for (let j = 0; j < Math.min(this._crt_bases.length, this.layers.length - 1); j++) {
      const layer = this.layers[j + 1]!;
      values = this._crt_bases[j]!.map((cs, i) => {
        let sum = 0n;
        for (let k = 0; k < cs.length && 2 * i + k < values.length; k++) {
          sum += cs[k]! * values[2 * i + k]!;
        }
        return remainder(sum, layer[i]!);
      });
    }
    if (values.length !== 1) throw new AssertionError('');
    return values[0]!;
  }
}

/**
 * Evaluate a product and its derivative with the original iterator_prod stack.
 * @see Reference: sage/rings/generic.py:prod_with_derivative; misc/misc_c.pyx:iterator_prod
 * @see Deviation: Product trees and factor-base nontermination
 */
export function prod_with_derivative(
  pairs: Iterable<readonly [IntegerLike, IntegerLike]>
): [bigint, bigint] {
  const combine = (a: [bigint, bigint], b: [bigint, bigint]): [bigint, bigint] => [
    a[0] * b[0],
    a[1] * b[0] + a[0] * b[1],
  ];
  const stack: Array<[bigint, bigint]> = [];
  let count = 0n;
  for (const [f, df] of pairs) {
    let value: [bigint, bigint] = [toBigInt(f), toBigInt(df)];
    count++;
    for (let n = count; (n & 1n) === 0n; n >>= 1n) value = combine(stack.pop()!, value);
    stack.push(value);
  }
  if (stack.length === 0)
    throw new TypeError("'sage.rings.integer.Integer' object is not iterable");
  while (stack.length > 1) {
    const right = stack.pop()!;
    stack.push(combine(stack.pop()!, right));
  }
  return stack[0]!;
}
