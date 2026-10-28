/**
 * Exact bounded small-scalar membership for the named-curve screening adapter.
 * Reuses cheatsheet arithmetic without exposing recovered scalars or pattern metadata.
 * @see Deviation: Weak-key application package (no SageMath counterpart)
 */
import {
  EcGroup,
  EcLanes,
  type EcPoint,
} from '../../zksecurity-cheatsheets/src/attacks/ec-arith.js';
import type { PrimeFieldCurveParams } from '../../zksecurity-cheatsheets/src/curves/params.js';

export function pointKey(point: EcPoint): string {
  return `${point.x.toString(16)}:${point.y.toString(16)}`;
}

/** Internal building block; the returned dictionary contains only public candidate points. */
export function buildSmallScalarDictionary(
  curve: PrimeFieldCurveParams,
  bound: number
): ReadonlySet<string> {
  if (!Number.isSafeInteger(bound) || bound < 1 || bound > 65536)
    throw new RangeError('Small-scalar bound must be between 1 and 65536');
  if (curve.order <= 2n * BigInt(bound))
    throw new RangeError('Small-scalar interval overlaps itself in this group');
  const group = new EcGroup({
    p: curve.field.p,
    a: curve.a,
    b: curve.b,
    n: curve.order,
    g: curve.generator,
  });
  if (!group.isOnCurve(group.g)) throw new RangeError('Generator is not on the curve');
  const block = new EcLanes([group.g]);
  while (block.length < bound) {
    const last = block.point(block.length - 1);
    if (last === null) throw new RangeError('Unexpected identity in candidate dictionary');
    const count = Math.min(block.length, bound - block.length);
    const next = new EcLanes([]);
    next.xs = block.xs.slice(0, count);
    next.ys = block.ys.slice(0, count);
    next.inf = block.inf.slice(0, count);
    next.add(group, new Array(count).fill(last.x), new Array(count).fill(last.y));
    if (next.inf.includes(1)) throw new RangeError('Unexpected identity in candidate dictionary');
    block.xs.push(...next.xs);
    block.ys.push(...next.ys);
    block.inf = new Uint8Array(block.xs.length);
  }
  const points = new Set<string>();
  for (let i = 0; i < bound; i++) {
    const point = block.point(i)!;
    if (!group.isOnCurve(point)) throw new RangeError('Candidate arithmetic verification failed');
    points.add(pointKey(point));
    points.add(pointKey(group.neg(point)!));
  }
  if (points.size !== 2 * bound) throw new RangeError('Candidate dictionary cardinality mismatch');
  return points;
}

// One bounded dictionary, not one allocation per key or every requested bound.
// Exact bound participates in the cache key: a prior wider search cannot widen a report.
let cached: { key: string; points: ReadonlySet<string> } | null = null;
export function hasSmallEcScalar(
  curve: PrimeFieldCurveParams,
  point: EcPoint,
  bound: number
): boolean {
  const key = [
    curve.field.p,
    curve.a,
    curve.b,
    curve.order,
    curve.generator.x,
    curve.generator.y,
    bound,
  ].join(':');
  if (cached?.key !== key) cached = { key, points: buildSmallScalarDictionary(curve, bound) };
  return cached.points.has(pointKey(point));
}
