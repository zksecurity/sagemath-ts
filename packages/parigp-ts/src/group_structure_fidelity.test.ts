import { expect, test } from 'bun:test';
import {
  group_elts,
  group_quotient,
  groupelts_quotient,
  quotient_perm,
  group_subgroups,
  cyclicgroup,
  trivialgroup,
} from './galconj.js';

// Native perm.c:962/1003 errors; six shared native comparisons cover callers.
test('quotient construction preserves the native unsupported-group error', () => {
  const G = cyclicgroup([0, 2, 1, 4, 3], 2);
  const H = cyclicgroup([0, 3, 4, 1, 2], 2);
  for (const run of [() => group_quotient(G, H), () => groupelts_quotient(group_elts(G, 4), H)]) {
    try {
      run();
      throw new Error('expected native rejection');
    } catch (e) {
      expect((e as Error).name).toBe('PariError');
      expect((e as Error).message).toBe(
        'sorry, group_quotient for a non-WSS group is not yet implemented.'
      );
    }
  }
});

test('quotient permutation preserves the native unsupported-group error', () => {
  const G = cyclicgroup([0, 2, 1, 4, 3], 2);
  const C = groupelts_quotient(group_elts(G, 4), trivialgroup());
  try {
    quotient_perm(C, [0, 3, 2, 1, 4]);
    throw new Error('expected native rejection');
  } catch (e) {
    expect((e as Error).name).toBe('PariError');
    expect((e as Error).message).toBe(
      'sorry, quotient_perm for a non-WSS group is not yet implemented.'
    );
  }
});

test('cyclic subgroups and quotient representatives match native order', () => {
  const G = cyclicgroup([0, 2, 3, 4, 1], 4);
  const H = cyclicgroup([0, 3, 4, 1, 2], 2);
  expect(group_subgroups(G)).toEqual([
    { gen: [[], H.gen[1]!, G.gen[1]!], ord: [0, 2, 2] },
    H,
    trivialgroup(),
  ]);
  expect(group_quotient(G, H)).toEqual({
    gen: [[], [0, 1, 2, 3, 4], [0, 2, 3, 4, 1]],
    coset: [0, 1, 2, 1, 2],
  });
});
