import { expect, test } from 'bun:test';
import { legacyNativeFixtures } from '../../../tests/property/native-live.mjs';
const fixtures = await legacyNativeFixtures(import.meta.url, './lll_real.fixtures.json');
import { fplll, fplll_heuristic } from './lll.js';
import { roundr_safe } from './gen3.js';
import { abscmprr, cmprr } from './kernel/none/cmp.js';
import type { MpReal } from './qfb.js';

type Fixture = {
  kind: string;
  probe: number;
  args: (number | string[])[];
  expected?: unknown;
  error?: string;
};
for (const fixture of fixtures as Fixture[]) {
  test(`native real LLL ${fixture.kind} probe ${fixture.probe}`, () => {
    const args = fixture.args,
      offset = fixture.kind === 'heuristic' ? 0 : 1;
    const m = Number(args[offset]),
      n = Number(args[offset + 1]),
      flat = args.at(-1) as string[];
    const original = Array.from({ length: n }, (_, j) =>
      Array.from({ length: m }, (_, i) => BigInt(flat[i * n + j]!))
    );
    const makeGram = (B: bigint[][]) =>
      B.map((a) => B.map((b) => a.reduce((sum, x, i) => sum + x * b[i]!, 0n)));
    const mode = offset ? Number(args[0]) : 0;
    const B = mode === 2 ? null : original.map((c) => [...c]);
    const G = mode ? makeGram(original) : null,
      savedG = G?.map((c) => [...c]) ?? null;
    let actual: unknown,
      error: unknown,
      U: bigint[][] | null = null,
      basis: bigint[][] | null = null,
      reducedG: bigint[][] | null = null;
    try {
      if (fixture.kind === 'heuristic') {
        const result = fplll_heuristic(
          B!,
          Number(args[2]) / Number(args[3]),
          Number(args[4]) / Number(args[5]),
          !!args[6],
          !!args[7],
          Number(args[8]),
          Number(args[9])
        );
        actual = result;
        [, basis, U] = result;
      } else {
        const [status, gram, reduced, transform, norms] = fplll(
          B,
          G,
          Number(args[3]) / Number(args[4]),
          Number(args[5]) / Number(args[6]),
          !!args[7],
          !!args[8],
          !!args[9],
          Number(args[10])
        );
        basis = reduced;
        reducedG = gram;
        U = transform;
        actual = [
          status,
          gram,
          reduced,
          transform,
          norms?.map((x) => [x.s, String(x.e), String(x.m), x.p]) ?? null,
        ];
      }
    } catch (caught) {
      error = caught;
    }
    expect(B).toEqual(mode === 2 ? null : original);
    expect(G).toEqual(savedG);
    if (fixture.error) {
      expect(error).toMatchObject({ name: 'PariError', message: fixture.error });
      return;
    }
    expect(error).toBeUndefined();
    expect(
      JSON.parse(JSON.stringify(actual, (_, x) => (typeof x === 'bigint' ? String(x) : x)))
    ).toEqual(fixture.expected);
    if (U) {
      const product = U.map((c) =>
        Array.from({ length: m }, (_, i) => c.reduce((sum, u, j) => sum + original[j]![i]! * u, 0n))
      );
      if (basis) expect(product).toEqual(basis);
      if (reducedG) {
        const exact = makeGram(product);
        for (let j = 0; j < n; j++)
          for (let i = 0; i <= j; i++) expect(reducedG[j]![i]).toBe(exact[j]![i]);
      }
    }
  });
}

test('roundr_safe preserves PARI half-integer ties and permits lost precision', () => {
  const x: MpReal = { s: 1, e: 1, m: 5n << 61n, p: 64 };
  expect(roundr_safe(x)).toBe(3n);
  expect(roundr_safe({ ...x, s: -1 })).toBe(-2n);
  expect(roundr_safe({ s: -1, e: -1, m: 1n << 63n, p: 64 })).toBe(0n);
  expect(roundr_safe({ s: 1, e: 128, m: 1n << 63n, p: 64 })).toBe(1n << 128n);
  expect(x).toEqual({ s: 1, e: 1, m: 5n << 61n, p: 64 });
});

test('absolute real comparison ignores zero accuracy and inspects trailing words', () => {
  const zero: MpReal = { s: 0, e: 1000, m: 0n, p: 0 },
    one: MpReal = { s: 1, e: 0, m: 1n << 63n, p: 64 };
  expect(cmprr(zero, one)).toBe(0);
  expect(abscmprr(zero, one)).toBe(-1);
  expect(abscmprr({ ...one, s: -1 }, one)).toBe(0);
  expect(abscmprr(one, { s: 1, e: 0, m: (1n << 127n) + 1n, p: 128 })).toBe(-1);
});

test('real LLL stage adapters enforce native precision and matrix preconditions', () => {
  for (const p of [0, 1, 63, 65, Infinity, NaN]) {
    expect(() => fplll_heuristic([[1n]], 0.99, 0.51, false, true, p)).toThrow(
      'LLL real precision must be a positive multiple of 64 bits'
    );
    expect(() => fplll_heuristic([[1n]], 0.99, 0.51, false, true, 64, p)).toThrow(
      'LLL real precision must be a positive multiple of 64 bits'
    );
    expect(() => fplll([[1n]], null, 0.99, 0.51, false, true, false, p)).toThrow(
      'LLL real precision must be a positive multiple of 64 bits'
    );
  }
  expect(() => fplll_heuristic([])).toThrow(
    'fplll_heuristic requires a nonempty rectangular basis'
  );
  expect(() => fplll(null)).toThrow(
    'fplll requires a nonempty rectangular basis or square Gram matrix'
  );
  expect(() =>
    fplll(
      [[1n]],
      [
        [1n, 0n],
        [0n, 1n],
      ]
    )
  ).toThrow('fplll requires a nonempty rectangular basis or square Gram matrix');
});
