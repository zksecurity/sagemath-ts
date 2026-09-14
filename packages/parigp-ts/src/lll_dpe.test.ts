import { expect, test } from 'bun:test';
import { legacyNativeFixtures } from '../../../tests/property/native-live.mjs';
const fixtures = await legacyNativeFixtures(import.meta.url, './lll_dpe.fixtures.json');
import { fplll_dpe } from './lll.js';

type Fixture = {
  probe: number;
  args: [number, number, number, number, number, number, number, number, number, number, string[]];
  expected?: unknown;
  error?: string;
  errorType?: string;
};

// Live native lll.c comparisons; no stored expected values.
for (const fixture of fixtures as unknown as Fixture[]) {
  test(`PARI DPE reduction state: native probe ${fixture.probe}`, () => {
    const [mode, m, n, dn, dd, en, ed, keep, track, want, flat] = fixture.args;
    const original = Array.from({ length: n }, (_, j) =>
      Array.from({ length: m }, (_, i) => BigInt(flat[i * n + j]!))
    );
    const gram = (B: bigint[][]) =>
      B.map((a) => B.map((b) => a.reduce((sum, x, i) => sum + x * b[i]!, 0n)));
    const B = mode === 2 ? null : original.map((c) => [...c]);
    const G = mode ? gram(original) : null;
    const savedG = G?.map((c) => [...c]) ?? null;
    let result: ReturnType<typeof fplll_dpe> | undefined;
    let error: unknown;
    try {
      result = fplll_dpe(B, G, dn / dd, en / ed, !!keep, !!track, !!want);
    } catch (caught) {
      error = caught;
    }
    expect(B).toEqual(mode === 2 ? null : original);
    expect(G).toEqual(savedG);
    if (fixture.error) {
      expect(error).toMatchObject({ name: fixture.errorType, message: fixture.error });
      return;
    }
    expect(error).toBeUndefined();
    const [status, reducedG, reducedB, U, norms] = result!;
    const state = [
      status,
      reducedG,
      reducedB,
      U,
      norms?.map((x) => [x.s, String(x.e), String(x.m), x.p]) ?? null,
    ];
    expect(
      JSON.parse(JSON.stringify(state, (_, v) => (typeof v === 'bigint' ? String(v) : v)))
    ).toEqual(fixture.expected);
    if (U) {
      const product = U.map((c) =>
        Array.from({ length: m }, (_, i) => c.reduce((sum, u, j) => sum + original[j]![i]! * u, 0n))
      );
      if (reducedB) expect(product).toEqual(reducedB);
      if (reducedG) {
        const exact = gram(product);
        // PARI updates only this triangle; the other half can retain old entries.
        for (let j = 0; j < n; j++)
          for (let i = 0; i <= j; i++) expect(reducedG[j]![i]).toBe(exact[j]![i]);
      }
    }
  });
}

test('DPE stage rejects adapters without a compatible native matrix shape', () => {
  for (const [B, G] of [
    [null, null],
    [[], null],
    [[[]], null],
    [[[1n], []], null],
    [null, [[1n], []]],
    [
      [[1n]],
      [
        [1n, 0n],
        [0n, 1n],
      ],
    ],
  ] as [bigint[][] | null, bigint[][] | null][]) {
    expect(() => fplll_dpe(B, G)).toThrow(
      'fplll_dpe requires a nonempty rectangular basis or square Gram matrix'
    );
  }
});
