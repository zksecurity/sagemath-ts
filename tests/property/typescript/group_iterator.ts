/** Fixtures for the bundled generic.py iterator and parent comparisons. */
import { ValueError } from '../../../packages/sagemath-ts/src/errors.js';
import { multiples, parseGroupOps } from '../../../packages/sagemath-ts/src/groups/generic.js';
import { QuadraticField } from '../../../packages/sagemath-ts/src/rings/number_field/number_field.js';

export function iterator_state(
  n: bigint,
  indexedFlag: bigint,
  copyKind: bigint,
  mutateBefore: bigint,
  mutateOp: bigint,
  throwAt: bigint
): string {
  const indexed = indexedFlag !== 0n;
  let copies = 0,
    calls = 0;
  class Box {
    declare __copy__?: () => Box | null;
    declare copy?: () => Box | null;
    value: number;
    constructor(value: number) {
      this.value = value % 11;
    }
    eq(other: Box): boolean {
      return this.value === other.value;
    }
  }
  if (copyKind !== 0n) {
    Box.prototype[copyKind === 4n ? 'copy' : '__copy__'] = function () {
      copies++;
      if (copyKind === 2n || (copyKind === 6n && this.value === 5))
        throw new ValueError('copy blocked');
      if (copyKind === 3n || (copyKind === 5n && this.value === 5)) return null;
      return new Box(this.value);
    };
  }
  const p = new Box(2),
    p0 = new Box(5),
    events: unknown[] = [];
  const state = () => [calls, copies, p.value, p0.value];
  const op = (x: Box, y: Box): Box => {
    calls++;
    if (BigInt(calls) === throwAt) throw new ValueError('step blocked');
    if (mutateOp !== 0n) {
      x.value = (x.value + y.value) % 11;
      return x;
    }
    return new Box(x.value + y.value);
  };
  let iterator: Generator<Box | [bigint, Box], void, undefined>;
  try {
    iterator = indexed
      ? multiples(p, n, p0, true, 'other', op)
      : multiples(p, n, p0, false, 'other', op);
    events.push(['construct', ...state()]);
  } catch (e) {
    events.push([
      'construct-error',
      e instanceof Error ? `${e.name}: ${e.message}` : String(e),
      ...state(),
    ]);
    return JSON.stringify(events);
  }
  if (mutateBefore !== 0n) {
    p.value = 9;
    p0.value = 8;
  }
  // The shared fixtures have n <= 3; this conversion only bounds observation.
  for (let i = 0; i < Math.max(Number(n), 0) + 2; i++) {
    try {
      const result = iterator.next();
      if (result.done) events.push(['done', ...state()]);
      else {
        const value = Array.isArray(result.value)
          ? [Number(result.value[0]), result.value[1].value]
          : result.value.value;
        events.push(['value', value, ...state()]);
      }
    } catch (e) {
      events.push([
        'next-error',
        e instanceof Error ? `${e.name}: ${e.message}` : String(e),
        ...state(),
      ]);
    }
  }
  return JSON.stringify(events);
}

export function iterator_standard(
  n: bigint,
  mode: bigint,
  supplied: bigint,
  parentKind: bigint,
  mutate: bigint
): string {
  const trace: unknown[] = [];
  type Parent = { zero(): Box; one(): Box };
  class Box {
    declare parent: Parent | (() => Parent);
    value: number;
    constructor(value: number) {
      this.value = value % 11;
    }
    eq(other: Box): boolean {
      return this.value === other.value;
    }
    add(other: Box): Box {
      trace.push(['add', this.value, other.value]);
      return new Box(this.value + other.value);
    }
    mul(other: Box | bigint): Box {
      trace.push([
        typeof other === 'bigint' ? 'scalar' : 'mul',
        this.value,
        typeof other === 'bigint' ? Number(other) : other.value,
      ]);
      return new Box(this.value * (typeof other === 'bigint' ? Number(other) : other.value));
    }
    pow(exponent: bigint): Box {
      trace.push(['pow', this.value, Number(exponent)]);
      return new Box(Number(BigInt(this.value) ** exponent % 11n));
    }
  }
  const ring: Parent = {
    zero() {
      trace.push(['zero']);
      return new Box(0);
    },
    one() {
      trace.push(['one']);
      return new Box(1);
    },
  };
  Box.prototype.parent =
    parentKind === 0n
      ? ring
      : function () {
          return ring;
        };
  const p = new Box(2),
    p0 = new Box(5);
  const iterator = multiples(
    p,
    n,
    supplied === 0n ? undefined : p0,
    false,
    mode === 0n ? '+' : '*'
  );
  const events: unknown[] = [['construct', trace.slice()]];
  if (mutate !== 0n) {
    p.value = 9;
    p0.value = 8;
  }
  for (let i = 0; i < Number(n) + 1; i++) {
    const value = iterator.next();
    events.push(
      value.done
        ? ['done', trace.slice(), p.value, p0.value]
        : ['value', value.value.value, trace.slice(), p.value, p0.value]
    );
  }
  return JSON.stringify(events);
}

export function parent_field(d: bigint, mode: bigint, k: bigint): string {
  const K = QuadraticField.create(d, 'a'),
    a = K.gen(),
    b = a.add(K.one());
  let result: unknown = null,
    error: string | null = null;
  try {
    const ops = parseGroupOps(mode === 0n ? '+' : '*', undefined, undefined, undefined, a);
    result = [
      ...[ops.identity, ops.inverse(a), ops.op(a, b), ops.power(a, k)].map((value) =>
        value.list().map(String)
      ),
      ops.isIdentity(ops.identity),
      ops.isIdentity(a),
    ];
  } catch (e) {
    error = e instanceof Error ? `${e.name}: ${e.message}` : String(e);
  }
  return JSON.stringify([result, error]);
}
