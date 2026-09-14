import { _nmod_poly_xgcd_kernels as k } from './gcd.js';
import { _nmod_poly_mul } from './mul.js';
import { _nmod_poly_evaluate_nmod } from './evaluate_nmod.js';
/** FLINT nmod_poly/compose.c and gr_poly composition/Taylor-shift dispatch.
 * @see Deviation: Polynomial Evaluation and Composition
 */
export function _nmod_poly_compose(
  a: readonly bigint[],
  b: readonly bigint[],
  p: bigint
): bigint[] {
  if (p < 2n) throw new RangeError('modulus must be at least 2');
  const A = k.normalized(a, p),
    B = k.normalized(b, p),
    norm = (a: readonly bigint[]) => k.normalized(a, p);
  const add = (a: readonly bigint[], b: readonly bigint[]) =>
    norm(
      Array.from({ length: Math.max(a.length, b.length) }, (_, i) => (a[i] ?? 0n) + (b[i] ?? 0n))
    );
  const mul = (a: readonly bigint[], b: readonly bigint[]) => _nmod_poly_mul(a, b, p);
  const horner = (a: readonly bigint[], b: readonly bigint[]) => {
    let out: bigint[] = [];
    for (let i = a.length - 1; i >= 0; i--) out = add(mul(out, b), [a[i]!]);
    return out;
  };
  const divide = (a: readonly bigint[], b: readonly bigint[]): bigint[] => {
    if (a.length <= 2 || b.length <= 1) return horner(a, b);
    let blocks: bigint[][] = [];
    for (let i = 0; i < a.length; i += 2)
      blocks.push(
        add(
          b.map((c) => c * (a[i + 1] ?? 0n)),
          [a[i]!]
        )
      );
    let power = mul(b, b);
    while (blocks.length > 2) {
      const next: bigint[][] = [];
      for (let i = 0; i < blocks.length; i += 2)
        next.push(i + 1 < blocks.length ? add(mul(power, blocks[i + 1]!), blocks[i]!) : blocks[i]!);
      blocks = next;
      power = mul(power, power);
    }
    return add(mul(power, blocks[1]!), blocks[0]!);
  };
  if (!A.length) return [];
  if (A.length === 1 || !B.length) return norm([A[0]!]);
  if (B.length === 1) return norm([_nmod_poly_evaluate_nmod(A, B[0]!, p)]);
  if (A.length <= 7) return horner(A, B);
  if (B.slice(1, -1).every((c) => c === 0n)) {
    let out = A.slice();
    const constant = B[0]!,
      scale = B[B.length - 1]!,
      degree = B.length - 1;
    if (constant !== 0n) {
      if (A.length <= 20) {
        for (let i = A.length - 2; i >= 0; i--)
          for (let j = i; j < A.length - 1; j++) out[j] = (out[j]! + constant * out[j + 1]!) % p;
      } else out = divide(A, [constant, 1n]);
    }
    if (scale === p - 1n) {
      for (let i = 1; i < out.length; i += 2) out[i] = out[i] === 0n ? 0n : p - out[i]!;
    } else if (scale !== 1n) {
      const bits = (out.length - 1).toString(2).length,
        powers = [scale];
      out[1] = (out[1]! * scale) % p;
      for (let j = 1; j < bits; j++) {
        powers.push((powers[j - 1]! * powers[j - 1]!) % p);
        out[2 ** j] = (out[2 ** j]! * powers[j]!) % p;
      }
      const zeros = (n: number) => {
        let c = 0;
        while (n % 2 === 0) {
          c++;
          n /= 2;
        }
        return c;
      };
      for (let i = 2 ** (bits - 1) + 1; i < out.length; i++) {
        const bit = zeros(i),
          slot = bits - 1 - bit;
        powers[slot] = (powers[slot]! * scale) % p;
        out[i / 2 ** bit] = (out[i / 2 ** bit]! * powers[slot]!) % p;
        for (let j = bit; j > 0; j--) {
          powers[bits - j] = (powers[bits - j - 1]! * powers[bits - j - 1]!) % p;
          out[i / 2 ** (j - 1)] = (out[i / 2 ** (j - 1)]! * powers[bits - j]!) % p;
        }
      }
      for (let i = Math.floor((out.length + 1) / 2); i < 2 ** (bits - 1); i++) {
        const bit = zeros(i),
          slot = bits - 2 - bit;
        powers[slot] = (powers[slot]! * scale) % p;
        out[i / 2 ** bit] = (out[i / 2 ** bit]! * powers[slot]!) % p;
        for (let j = bit; j > 0; j--) {
          powers[bits - j - 1] = (powers[bits - j - 2]! * powers[bits - j - 2]!) % p;
          out[i / 2 ** (j - 1)] = (out[i / 2 ** (j - 1)]! * powers[bits - j - 1]!) % p;
        }
      }
    }
    if (degree > 1) {
      const stretched = Array<bigint>((out.length - 1) * degree + 1).fill(0n);
      for (let i = 0; i < out.length; i++) stretched[i * degree] = out[i]!;
      out = stretched;
    }
    return norm(out);
  }
  return divide(A, B);
}
