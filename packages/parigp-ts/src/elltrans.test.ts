import { expect, test } from 'bun:test';
import { sumdedekind } from './elltrans.js';

// Observed through cypari2; shared comparative cases permanently repeat these calls.
const cases = [
  {
    function: 'dedekind_native',
    args: ['-5', '-6148914691236517205', '1'],
    result: '[126030506267014245406065791099044180, 1229782938247303441]',
    error: null,
    errorType: null,
    seed: 1,
  },
  {
    function: 'dedekind_native',
    args: ['1', '-6148914691236517205', '1'],
    result: '[-6301525313350712288750033628661760607, 12297829382473034410]',
    error: null,
    errorType: null,
    seed: 1,
  },
  {
    function: 'dedekind_native',
    args: ['5', '-6148914691236517205', '1'],
    result: '[-126030506267014245406065791099044180, 1229782938247303441]',
    error: null,
    errorType: null,
    seed: 1,
  },
  {
    function: 'dedekind_native',
    args: ['-5', '-6148914691236517204', '1'],
    result: '[1260305062670142456110296141402614201, 12297829382473034408]',
    error: null,
    errorType: null,
    seed: 1,
  },
  {
    function: 'dedekind_native',
    args: ['1', '-6148914691236517204', '1'],
    result: '[-6301525313350712280551480707013071001, 12297829382473034408]',
    error: null,
    errorType: null,
    seed: 1,
  },
  {
    function: 'dedekind_native',
    args: ['5', '-6148914691236517204', '1'],
    result: '[-1260305062670142456110296141402614201, 12297829382473034408]',
    error: null,
    errorType: null,
    seed: 1,
  },
  {
    function: 'dedekind_native',
    args: ['-5', '0', '1'],
    result: null,
    error: 'impossible inverse in dvmdii: 0',
    errorType: 'PariError',
    seed: 1,
  },
  {
    function: 'dedekind_native',
    args: ['1', '0', '1'],
    result: null,
    error: 'impossible inverse in dvmdii: 0',
    errorType: 'PariError',
    seed: 1,
  },
  {
    function: 'dedekind_native',
    args: ['5', '0', '1'],
    result: null,
    error: 'impossible inverse in dvmdii: 0',
    errorType: 'PariError',
    seed: 1,
  },
  {
    function: 'dedekind_native',
    args: ['-5', '9223372036854775808', '1'],
    result: '[-2835686391007820524249769043504013859, 18446744073709551616]',
    error: null,
    errorType: null,
    seed: 1,
  },
  {
    function: 'dedekind_native',
    args: ['1', '9223372036854775808', '1'],
    result: '[14178431955039102639695589291229620907, 18446744073709551616]',
    error: null,
    errorType: null,
    seed: 1,
  },
  {
    function: 'dedekind_native',
    args: ['5', '9223372036854775808', '1'],
    result: '[2835686391007820524249769043504013859, 18446744073709551616]',
    error: null,
    errorType: null,
    seed: 1,
  },
];
for (const row of cases)
  test(`PARI original Dedekind boundary ${row.args[0]}/${row.args[1]}`, () => {
    const call = () => sumdedekind(BigInt(row.args[0]!), BigInt(row.args[1]!));
    if (row.error !== null) {
      let error: Error | undefined;
      try {
        call();
      } catch (e) {
        error = e as Error;
      }
      expect(error?.name).toBe(row.errorType);
      expect(error?.message).toBe(row.error);
    } else {
      const result = call();
      expect(`[${result.join(', ')}]`).toBe(row.result!);
    }
  });
