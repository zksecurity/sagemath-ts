/** Real-double scalar coercion and generic multiple operation schedules. */
import { RDF } from '../../../packages/sagemath-ts/src/rings/real_double.js';
import { Integer } from '../../../packages/sagemath-ts/src/rings/integer_ring.js';
import { Rational } from '../../../packages/sagemath-ts/src/rings/rational.js';
import { multiple } from '../../../packages/sagemath-ts/src/groups/generic.js';
function fromBits(bits: bigint): number {
  const v = new DataView(new ArrayBuffer(8));
  v.setBigUint64(0, bits);
  return v.getFloat64(0);
}
function bits(x: number): string {
  if (Number.isNaN(x)) return 'NaN';
  const v = new DataView(new ArrayBuffer(8));
  v.setFloat64(0, x);
  return String(v.getBigUint64(0));
}
export function rdf_scalar(
  fn: bigint,
  left: bigint,
  n: bigint,
  d: bigint,
  right: bigint,
  kind: bigint
): string {
  const a = RDF.__call__(fromBits(left));
  const b =
    kind === 0n
      ? n
      : kind === 1n
        ? new Integer(n)
        : kind === 2n
          ? new Rational(n, d)
          : kind === 3n
            ? fromBits(right)
            : kind === 4n
              ? RDF.__call__(fromBits(right))
              : n !== 0n;
  const methods = [() => a.add(b), () => a.sub(b), () => a.mul(b), () => a.div(b), () => a.eq(b)];
  const r = methods[Number(fn)]!();
  return typeof r === 'boolean' ? String(r) : bits(r.value);
}
export function rdf_multiple(value: bigint, n: bigint, mode: bigint): string {
  return bits(multiple(RDF.__call__(fromBits(value)), n, mode === 0n ? '+' : '*').value);
}

export function rdf_alias(fn: bigint, left: bigint, n: bigint, kind: bigint, field = RDF): string {
  const a = field.__call__(fromBits(left));
  const b =
    kind === 0n
      ? n
      : kind === 1n
        ? new Integer(n)
        : kind === 2n
          ? new Rational(n)
          : kind === 3n
            ? Number(n)
            : kind === 4n
              ? field.__call__(n)
              : n !== 0n;
  const methods = [() => a.add(b), () => a.sub(b), () => a.mul(b), () => a.div(b)];
  const r = methods[Number(fn)]!();
  return JSON.stringify([bits(r.value), r === a]);
}
