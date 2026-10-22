import { describe, expect, test } from 'bun:test';
import { createHash } from 'node:crypto';
import { sqrt_mod } from '../../../sagemath-ts/src/arith/misc.js';
import {
  type PrimeFieldCurveParams,
  bls12_381_g1,
  bn254_g1,
  createPrimeFieldCurve,
  secp256k1,
} from '../curves.js';
import { EcGroup, EcLanes, type EcPoint, ecPlusMinusX } from './ec-arith.js';
import {
  ecCurveReport,
  findCloseEcKeys,
  findEcKeysInInterval,
  findPatternedEcKeys,
  findSmallEcKeys,
  findWeakEcKeys,
  validateEcKeys,
} from './ec-keys.js';
import { EcBabyTable } from './ec-table.js';

const p384: PrimeFieldCurveParams = {
  name: 'secp384r1',
  field: {
    p: 0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffeffffffff0000000000000000ffffffffn,
  },
  a: -3n,
  b: 0xb3312fa7e23ee7e4988e056be3f82d19181d9c6efe8141120314088f5013875ac656398d8a2ed19d2a85c8edd3ec2aefn,
  order:
    0xffffffffffffffffffffffffffffffffffffffffffffffffc7634d81f4372ddf581a0db248b0a77aecec196accc52973n,
  cofactor: 1n,
  generator: {
    x: 0xaa87ca22be8b05378eb1c71ef320ad746e1d3b628ba79b9859f741e082542a385502f25dbf55296c3a545e3872760ab7n,
    y: 0x3617de4a96262c6f5d9e98bf9292dc29f8f41dbd289a147ce9da3113b5f0b8c00a60b1ce1d7e819d7a431d7c90ea0e5fn,
  },
};

function groupOf(curve: PrimeFieldCurveParams): EcGroup {
  return new EcGroup({
    p: curve.field.p,
    a: curve.a,
    b: curve.b,
    n: curve.order,
    g: curve.generator,
  });
}

/** Deterministic "sound" scalars: SHA-256 chains reduced mod n. */
function randomScalar(label: string, n: bigint): bigint {
  const h1 = createHash('sha256').update(`${label}/1`).digest('hex');
  const h2 = createHash('sha256').update(`${label}/2`).digest('hex');
  return BigInt(`0x${h1}${h2}`) % n;
}

function pointAt(p: bigint, a: bigint, b: bigint, from: bigint): EcPoint {
  for (let x = from; ; x++) {
    const y = sqrt_mod((x * x * x + a * x + b) % p, p);
    if (y !== null) return { x, y };
  }
}

function repeatedWord(word: bigint, words: number): bigint {
  let s = 0n;
  for (let i = 0; i < words; i++) s = (s << 32n) | word;
  return s;
}

describe('EcGroup arithmetic', () => {
  test('scalar multiplication matches sagemath-ts on secp256k1', () => {
    const group = groupOf(secp256k1);
    const { generator } = createPrimeFieldCurve(secp256k1);
    for (const k of [1n, 2n, 3n, 0xdeadbeefn, randomScalar('mul', secp256k1.order)]) {
      const ours = group.mul(group.g, k) as EcPoint;
      const theirs = generator.mul(k);
      expect(ours).toEqual({ x: BigInt(String(theirs.x)), y: BigInt(String(theirs.y)) });
    }
    expect(group.mul(group.g, secp256k1.order)).toBeNull();
    expect(group.mul(group.g, -5n)).toEqual(group.neg(group.mul(group.g, 5n)));
  });

  test('batched lanes agree with single additions, including doubling and cancellation', () => {
    const group = groupOf(p384);
    const P = group.mul(group.g, 7n) as EcPoint;
    const Q = group.mul(group.g, 11n) as EcPoint;
    const lanes = new EcLanes([P, Q, group.neg(Q), null, group.mul(group.g, 99n)]);
    lanes.add(group, Array(5).fill(Q.x), Array(5).fill(Q.y));
    expect(lanes.point(0)).toEqual(group.mul(group.g, 18n));
    expect(lanes.point(1)).toEqual(group.mul(group.g, 22n));
    expect(lanes.point(2)).toBeNull();
    expect(lanes.point(3)).toEqual(Q);
    expect(lanes.point(4)).toEqual(group.mul(group.g, 110n));

    const blockPoints = [1n, 2n, 3n, 9n].map((k) => group.mul(group.g, k) as EcPoint);
    const plus: bigint[] = [];
    const minus: bigint[] = [];
    ecPlusMinusX(
      group,
      new EcLanes([group.mul(group.g, 9n), null]),
      blockPoints.map((B) => B.x),
      blockPoints.map((B) => B.y),
      plus,
      minus
    );
    expect(plus.slice(0, 3)).toEqual(
      [10n, 11n, 12n].map((k) => (group.mul(group.g, k) as EcPoint).x)
    );
    expect(minus.slice(0, 3)).toEqual(
      [8n, 7n, 6n].map((k) => (group.mul(group.g, k) as EcPoint).x)
    );
    expect(plus[3]).toBe((group.mul(group.g, 18n) as EcPoint).x);
    expect(minus[3]).toBe(-1n);
    expect(plus.slice(4)).toEqual(blockPoints.map((B) => B.x));
  });

  test('the baby-step table answers ±j for every stored multiple', () => {
    const group = groupOf(bn254_g1);
    const table = new EcBabyTable(group, group.g, 12);
    for (const j of [1, 2, 2048, 2049, 4095, 4096]) {
      expect(table.get((group.mul(group.g, BigInt(j)) as EcPoint).x)).toBe(j);
      expect(table.get((group.mul(group.g, -BigInt(j)) as EcPoint).x)).toBe(j);
    }
    expect(table.get((group.mul(group.g, 4097n) as EcPoint).x)).toBe(0);
    expect(table.stride).toEqual(group.mul(group.g, 8193n) as EcPoint);
  });
});

describe('ecCurveReport', () => {
  test('standard curves have no problems or weaknesses', () => {
    for (const curve of [secp256k1, p384, bn254_g1]) {
      const report = ecCurveReport(curve);
      expect(report.problems).toEqual([]);
      expect(report.weaknesses).toEqual([]);
      expect(report.largestPrimeFactor).toBe(curve.order);
      expect(report.anomalous).toBe(false);
    }
    expect(ecCurveReport(secp256k1).embeddingDegree).toBeNull();
    expect(ecCurveReport(bn254_g1).embeddingDegree).toBe(12);
  });

  test('the Hasse bound forces the BLS12-381 G1 cofactor', () => {
    const report = ecCurveReport(bls12_381_g1);
    expect(report.cofactor).toBe(0x396c8c005555e1568c00aaab0000aaabn);
    expect(report.embeddingDegree).toBe(12);
    expect(report.weaknesses).toEqual([]);
  });

  test('a supersingular curve has embedding degree 2', () => {
    const ell = 1208925819614629174706261n;
    const p = 4n * ell - 1n;
    const group = new EcGroup({ p, a: 1n, b: 0n, n: ell, g: { x: 0n, y: 0n } });
    const generator = group.mul(pointAt(p, 1n, 0n, 2n), 4n) as EcPoint;
    const report = ecCurveReport({ name: 'ss', field: { p }, a: 1n, b: 0n, order: ell, generator });
    expect(report.problems).toEqual([]);
    expect(report.embeddingDegree).toBe(2);
    expect(report.curveOrder).toBe(p + 1n);
    expect(report.weaknesses.join()).toContain('embedding degree 2');
  });

  test("an anomalous curve is flagged for Smart's attack", () => {
    const p = 730750818665451459112596905638433048232067471723n;
    const a = 425706413842211054102700238164133538302169176474n;
    const b = 203362936548826936673264444982866339953265530166n;
    const report = ecCurveReport({
      name: 'anomalous',
      field: { p },
      a,
      b,
      order: p,
      generator: pointAt(p, a, b, 2n),
    });
    expect(report.problems).toEqual([]);
    expect(report.anomalous).toBe(true);
    expect(report.curveOrder).toBe(p);
  });

  test('a small or composite order is reported with its factors', () => {
    const p = 18446744073709551557n;
    const n = 18446744070893002621n;
    const generator = pointAt(p, 3n, 10n, 1n);
    const small = ecCurveReport({ name: 'toy', field: { p }, a: 3n, b: 10n, order: n, generator });
    expect(small.problems).toEqual([]);
    expect(small.weaknesses.join()).toContain('64 bits');

    const doubled = ecCurveReport({
      name: 'bad',
      field: { p },
      a: 3n,
      b: 10n,
      order: 6n * n,
      generator,
    });
    expect(doubled.problems).toContain('order n is not prime; pass the prime subgroup order');
    expect(doubled.orderFactors).toEqual([
      [2n, 1n],
      [3n, 1n],
      [n, 1n],
    ]);
  });

  test('a singular curve and a wrong generator are caught', () => {
    const p = 1000003n;
    const cusp = ecCurveReport({
      name: 'cusp',
      field: { p },
      a: 0n,
      b: 0n,
      order: p,
      generator: { x: 1n, y: 1n },
    });
    expect(cusp.weaknesses.some((w) => w.startsWith('singular'))).toBe(true);
    expect(cusp.anomalous).toBe(true);
    const off = ecCurveReport({ ...secp256k1, generator: { x: 1n, y: 1n } });
    expect(off.problems).toContain('generator is not on the curve');
  });
});

describe('validateEcKeys', () => {
  test('decodes SEC1 encodings and rejects bad keys', () => {
    const group = groupOf(secp256k1);
    const Q = group.mul(group.g, 12345n) as EcPoint;
    const hex = (v: bigint) => v.toString(16).padStart(64, '0');
    const compressed = `0${2n + (Q.y & 1n)}${hex(Q.x)}`;
    const checks = validateEcKeys(secp256k1, [
      Q,
      `04${hex(Q.x)}${hex(Q.y)}`,
      compressed,
      Uint8Array.from(Buffer.from(compressed, 'hex')),
      { x: Q.x, y: Q.y + 1n },
      '00',
      'zz',
      `02${'ff'.repeat(32)}`,
    ]);
    expect(checks.slice(0, 4).map((c) => c.point)).toEqual([Q, Q, Q, Q]);
    expect(checks.map((c) => c.problem)).toEqual([
      undefined,
      undefined,
      undefined,
      undefined,
      'not-on-curve',
      'identity',
      'malformed',
      'not-on-curve',
    ]);
  });

  test('detects points outside the prime-order subgroup', () => {
    const report = ecCurveReport(bls12_381_g1);
    const outside = pointAt(bls12_381_g1.field.p, 0n, 4n, 5n);
    const inside = groupOf(bls12_381_g1).mul(outside, report.cofactor as bigint) as EcPoint;
    const checks = validateEcKeys(bls12_381_g1, [outside, inside]);
    expect(checks.map((c) => c.problem)).toEqual(['wrong-subgroup', undefined]);
  });
});

describe('key searches', () => {
  test('patterned P-384 scalars are recovered and sound keys are not', () => {
    const n = p384.order;
    const group = groupOf(p384);
    const secrets = [
      255n,
      randomScalar('p384-a', n),
      (1n << 20n) + 7n,
      0xdeadbeefn,
      n - 5n,
      0xabcdefn << 160n,
      repeatedWord(0xcafebaben, 12),
      repeatedWord(0x01020304n, 3),
      randomScalar('p384-b', n),
    ];
    const keys = secrets.map((s) => group.mul(group.g, s) as EcPoint);
    const found = findPatternedEcKeys(p384, keys);
    expect(found.map((f) => f.index)).toEqual([0, 2, 3, 4, 5, 6, 7]);
    for (const f of found) expect(f.scalar).toBe(secrets[f.index] as bigint);
    expect(found.map((f) => f.pattern)).toEqual([
      { kind: 'small', value: 255n, negated: false },
      { kind: 'small', value: (1n << 20n) + 7n, negated: false },
      { kind: 'small', value: 0xdeadbeefn, negated: false },
      { kind: 'small', value: 5n, negated: true },
      { kind: 'shifted', value: 0xabcdefn, shift: 20, negated: false },
      { kind: 'repeated-word', word: 0xcafebaben, words: 12, negated: false },
      { kind: 'repeated-word', word: 0x01020304n, words: 3, negated: false },
    ]);

    // Coverage grows with whatever table is cached; it always includes ±2^16.
    const small = findSmallEcKeys(p384, keys).map((f) => f.index);
    expect(small).toContain(0);
    expect(small).toContain(4);
    expect(small).not.toContain(1);
    expect(small).not.toContain(8);
  });

  test('close pairs and repeated keys on BN254', () => {
    const n = bn254_g1.order;
    const group = groupOf(bn254_g1);
    const s = randomScalar('bn-close', n);
    const secrets = [s, s + 1000n, s - 9_000_000n, s, randomScalar('bn-far', n)];
    const keys = secrets.map((k) => group.mul(group.g, k) as EcPoint);
    expect(findCloseEcKeys(bn254_g1, keys)).toEqual([
      { a: 0, b: 1, difference: -1000n },
      { a: 0, b: 2, difference: 9_000_000n },
      { a: 0, b: 3, difference: 0n },
      { a: 1, b: 2, difference: 9_001_000n },
      { a: 1, b: 3, difference: 1000n },
      { a: 2, b: 3, difference: -9_000_000n },
    ]);
  });

  test('the kangaroo finds scalars in an interval', () => {
    const n = bn254_g1.order;
    const group = groupOf(bn254_g1);
    const secrets = [0x2_9876_5432n, randomScalar('bn-kangaroo', n), 0x3_0000_0001n];
    const keys = secrets.map((k) => group.mul(group.g, k) as EcPoint);
    const result = findEcKeysInInterval(bn254_g1, keys, { bits: 34 });
    expect(result.recovered).toEqual([
      { index: 0, scalar: secrets[0] as bigint, pattern: { kind: 'interval', bits: 34 } },
      { index: 2, scalar: secrets[2] as bigint, pattern: { kind: 'interval', bits: 34 } },
    ]);
    expect(result.steps).toBeLessThanOrEqual(Math.ceil(6 * Math.sqrt(3 * 2 ** 34)) + 1024);
  });

  test('findWeakEcKeys combines the checks and derives related keys', () => {
    const n = secp256k1.order;
    const group = groupOf(secp256k1);
    const shifted = 0x1234n << 200n;
    const secrets = [shifted, shifted + 77_777n, randomScalar('k1', n), 3n];
    const keys: (EcPoint | string)[] = secrets.map((k) => group.mul(group.g, k) as EcPoint);
    keys.push('00', { x: 1n, y: 1n } as EcPoint);
    const report = findWeakEcKeys(secp256k1, keys);
    expect(report.curve.weaknesses).toEqual([]);
    expect(report.invalid.map((c) => [c.index, c.problem])).toEqual([
      [4, 'identity'],
      [5, 'not-on-curve'],
    ]);
    expect(
      report.recovered.map((r) => [r.index, r.scalar === secrets[r.index], r.pattern.kind])
    ).toEqual([
      [0, true, 'shifted'],
      [1, true, 'related'],
      [3, true, 'small'],
    ]);
    expect(report.close).toEqual([{ a: 0, b: 1, difference: -77_777n }]);
  });
});
