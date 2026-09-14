import { expect, test } from 'bun:test';
import { listznstarelts } from './galconj.js';
import { znstar } from './char.js';
import { subgrouplist } from './subgroup.js';
import { ZM_hnfmodid } from './hnf_snf.js';
import { pgener_Zp } from './arith1.js';

test('unit subgroups retain native Birkhoff ordering', () => {
  expect(listznstarelts(8, 4)).toEqual([
    [0, 1],
    [0, 1, 3],
    [0, 1, 7],
    [0, 1, 5],
    [0, 1, 3, 5, 7],
  ]);
  expect(listznstarelts(8, 2)).toEqual([
    [0, 1],
    [0, 1, 3],
    [0, 1, 7],
    [0, 1, 5],
  ]);
});

test('cyclic decomposition uses native Conrey generators', () => {
  expect(znstar(8n)).toEqual([4n, [2n, 2n], [7n, 5n]]);
  expect(znstar(0n)).toEqual([2n, [2n], [-1n]]);
  expect(pgener_Zp(40487n)).toBe(10n);
});

test('exact-index enumeration retains native multiword packing', () => {
  expect(subgrouplist([1n << 65n], [1n << 65n])).toEqual([[[2n]]]);
  expect(subgrouplist([8n], [2n])).toEqual([[[2n]]]);
  expect(subgrouplist([8n], [3n])).toEqual([]);
});

test('modular HNF inserts missing diagonal columns and reduces signed entries', () => {
  expect(ZM_hnfmodid([], [6n, 4n])).toEqual([
    [6n, 0n],
    [0n, 4n],
  ]);
  expect(ZM_hnfmodid([[-2n, 0n]], [6n, 4n])).toEqual([
    [2n, 0n],
    [0n, 4n],
  ]);
});

test('new dependencies preserve native domain and dimension errors', () => {
  const cases: [() => unknown, string][] = [
    [() => pgener_Zp(2n), 'domain error in pgener_Zl: p = 2'],
    [
      () => subgrouplist([2n, 3n], [1n]),
      'incorrect type in subgrouplist [not a finite group] (t_VEC).',
    ],
    [() => subgrouplist([2n], [0n]), 'domain error in subgroup: index bound <= 0'],
    [() => ZM_hnfmodid([[1n, 2n]], [2n]), 'inconsistent dimensions in ZM_hnfmod.'],
  ];
  for (const [run, message] of cases) {
    try {
      run();
      throw new Error('expected native rejection');
    } catch (e) {
      expect((e as Error).name).toBe('PariError');
      expect((e as Error).message).toBe(message);
    }
  }
});

test('negative order arguments use native unsigned-word conversion', () => {
  expect(listznstarelts(9, -1)).toEqual([
    [0, 1],
    [0, 1, 4, 7],
  ]);
  expect(listznstarelts(9, -3)).toEqual([[0, 1]]);
});
