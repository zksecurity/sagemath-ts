import { generic_power } from '../../../packages/sagemath-ts/src/arith/power.js';
import { Integer } from '../../../packages/sagemath-ts/src/rings/integer_ring.js';
import { ZeroDivisionError } from '../../../packages/sagemath-ts/src/errors.js';
export function power_case(value: bigint, exponent: bigint, kind: bigint): string {
  const trace: unknown[][] = [],
    ring = {
      one: () => {
        trace.push(['one']);
        return new Box(1n);
      },
    };
  class Box {
    v: bigint;
    constructor(v: bigint) {
      this.v = ((v % 101n) + 101n) % 101n;
    }
    parent() {
      return ring;
    }
    mul(b: Box) {
      trace.push(['mul', Number(this.v), Number(b.v)]);
      return new Box(this.v * b.v);
    }
    inv() {
      trace.push(['inv', Number(this.v)]);
      if (this.v === 0n) throw new ZeroDivisionError('zero inverse');
      let r = 1n,
        b = this.v,
        n = 99n;
      while (n) {
        if (n & 1n) r = (r * b) % 101n;
        b = (b * b) % 101n;
        n >>= 1n;
      }
      return new Box(r);
    }
  }
  const a = new Box(value);
  let result = null,
    error = null;
  try {
    const r = generic_power(a, kind === 0n ? exponent : new Integer(exponent));
    result = [String(r.v), r === a];
  } catch (e) {
    error = e.name + ': ' + e.message;
  }
  return JSON.stringify([result, error, trace]);
}
