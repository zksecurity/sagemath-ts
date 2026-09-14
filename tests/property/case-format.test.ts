import { describe, expect, test } from 'bun:test';
import { readdirSync } from 'node:fs';
import { join } from 'node:path';
import {
  type CaseSuite,
  countCases,
  decodeRow,
  encodeScalar,
  isCanonicalText,
  isNormalizedShape,
  normalizeSuite,
  parseFixedGenerator,
  serializeSuiteToString,
} from './case-format.js';

import { isCaseFile, readText } from './storage.js';

const CASES_DIR = join(import.meta.dir, 'cases');

describe('case-format scalars', () => {
  test('safe integers round-trip as numbers, larger ones as strings', () => {
    expect(encodeScalar(5n)).toBe(5);
    expect(encodeScalar(-(2n ** 53n) + 1n)).toBe(-(2 ** 53) + 1);
    expect(encodeScalar(2n ** 53n)).toBe('9007199254740992');
    expect(decodeRow([1, 5, '9007199254740992', [1, '-9007199254740992']])).toEqual({
      seed: 1,
      args: [5n, 9007199254740992n, [1n, -9007199254740992n]],
    });
  });

  test('unsafe JSON numbers and non-integer text are rejected', () => {
    expect(() => decodeRow([1, 2 ** 53])).toThrow('Unsafe integer');
    expect(() => decodeRow([1, '1.5'])).toThrow('Invalid integer');
    expect(() => decodeRow([1.5, 1])).toThrow('Invalid seed');
  });

  test('parseFixedGenerator handles ints, lists, empty lists, and literals', () => {
    expect(parseFixedGenerator('fixedValue(7)')).toBe(7n);
    expect(parseFixedGenerator('fixedValue([1, -2,3])')).toEqual([1n, -2n, 3n]);
    expect(parseFixedGenerator('fixedValue([])')).toEqual([]);
    expect(parseFixedGenerator('-12')).toBe(-12n);
    expect(parseFixedGenerator('randomBigint(1, 10)')).toBeUndefined();
    expect(parseFixedGenerator('fixedValue(abc)')).toBeUndefined();
  });
});

describe('normalizeSuite', () => {
  const verbose: CaseSuite = {
    module: 'demo',
    description: 'd',
    cases: [
      {
        function: 'f',
        seeds: [1],
        argGenerators: ['fixedValue(11)', 'fixedValue([1, 1])'],
        description: 'x',
      },
      { function: 'f', seeds: [2, 3], argGenerators: ['fixedValue(7)', 'fixedValue([])'] },
      { function: 'f', seeds: [1], argGenerators: ['fixedValue(11)', 'fixedValue([1, 1])'] },
      { function: 'g', seeds: [4, 5], argGenerators: ['randomBigint(1, 10)'], _comment: 'keep' },
      { function: 'g', seeds: [4, 5], argGenerators: ['randomBigint(1, 10)'], _comment: 'keep' },
      { function: 'f', rows: [[9, '9007199254740993', [2]]] },
      { function: 'f', rows: [[9, '9007199254740993', [2]]] },
    ],
  };

  test('converts fixed cases to rows, merges neighbours, dedupes, keeps generators', () => {
    const suite = normalizeSuite(verbose);
    expect(suite).toEqual({
      module: 'demo',
      description: 'd',
      cases: [
        {
          function: 'f',
          rows: [
            [1, 11, [1, 1]],
            [2, 7, []],
            [3, 7, []],
          ],
        },
        { function: 'g', seeds: [4, 5], argGenerators: ['randomBigint(1, 10)'], _comment: 'keep' },
        { function: 'f', rows: [[9, '9007199254740993', [2]]] },
      ],
    });
    expect(countCases(suite)).toBe(6);
    expect(countCases(verbose)).toBe(10);
  });

  test('serialization is stable and canonical', () => {
    const text = serializeSuiteToString(normalizeSuite(verbose));
    expect(text).toBe(
      [
        '{',
        '  "module": "demo",',
        '  "description": "d",',
        '  "cases": [',
        '    {',
        '      "function": "f",',
        '      "rows": [',
        '        [1, 11, [1, 1]],',
        '        [2, 7, []],',
        '        [3, 7, []]',
        '      ]',
        '    },',
        '    {',
        '      "function": "g",',
        '      "seeds": [',
        '        4,',
        '        5',
        '      ],',
        '      "argGenerators": [',
        '        "randomBigint(1, 10)"',
        '      ],',
        '      "_comment": "keep"',
        '    },',
        '    {',
        '      "function": "f",',
        '      "rows": [',
        '        [9, "9007199254740993", [2]]',
        '      ]',
        '    }',
        '  ]',
        '}',
        '',
      ].join('\n')
    );
    expect(JSON.parse(text)).toEqual(normalizeSuite(verbose));
    expect(isCanonicalText(text)).toBe(true);
    expect(isCanonicalText(JSON.stringify(verbose))).toBe(false);
    expect(isCanonicalText(serializeSuiteToString(verbose))).toBe(false);
  });

  test('isNormalizedShape rejects what normalizeSuite would change', () => {
    const ok = (cases: CaseSuite['cases']) => isNormalizedShape({ module: 'm', cases });
    expect(ok([{ function: 'f', rows: [[1, 2, ['9007199254740992']]] }])).toBe(true);
    expect(ok([{ function: 'f', rows: [[1, '5']] }])).toBe(false);
    expect(ok([{ function: 'f', rows: [[1, '9007199254740991']] }])).toBe(false);
    expect(ok([{ function: 'f', rows: [[1, 2 ** 53]] }])).toBe(false);
    expect(ok([{ function: 'f', rows: [] }])).toBe(false);
    expect(
      ok([
        { function: 'f', rows: [[1]] },
        { function: 'f', rows: [[2]] },
      ])
    ).toBe(false);
    expect(ok([{ function: 'f', seeds: [1], argGenerators: ['fixedValue(1)'] }])).toBe(false);
    expect(ok([{ function: 'f', seeds: [1], argGenerators: ['randomBigint(1, 9)'] }])).toBe(true);
  });

  test('empty suite serializes', () => {
    const text = serializeSuiteToString({ module: 'e', cases: [] });
    expect(text).toBe('{\n  "module": "e",\n  "cases": []\n}\n');
    expect(isCanonicalText(text)).toBe(true);
  });
});

describe('checked-in case files', () => {
  const files = readdirSync(CASES_DIR)
    .filter(isCaseFile)
    .sort();

  test.each(files)(
    '%s is canonical (else run `bun tests/property/normalize-cases.ts`)',
    (file) => {
      const text = readText(join(CASES_DIR, file));
      expect(isCanonicalText(text)).toBe(true);
    },
    60_000
  );
});
