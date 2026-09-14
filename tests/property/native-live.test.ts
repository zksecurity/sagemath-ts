import { expect, test } from 'bun:test';
import { readFileSync } from 'node:fs';
import { nativeInputs } from './native-live.mjs';

const catalog = JSON.parse(readFileSync(new URL('./native-suites.json', import.meta.url), 'utf8'));
const names = Object.keys(catalog);
const payloads = (name: string, seed: number) => nativeInputs(name, seed, 8)
  .flatMap(s => s.cases.flatMap(c => c.rows.map(r => r.slice(1))));

test('every native profile materializes compact, reproducible inputs without saved corpora', () => {
  expect(names).toHaveLength(97);
  let regressions = 0;
  for (const name of names) {
    const suites = nativeInputs(name, 123, 8);
    expect(suites, name).toEqual(nativeInputs(name, 123, 8));
    for (const suite of suites) {
      expect(suite.cases[0]!.rows.length).toBeGreaterThan(0);
      if (suite.regressionId !== undefined) {
        regressions++;
        expect(suite).toEqual(nativeInputs(name, 456, 8).find(s => s.regressionId === suite.regressionId));
      }
    }
    expect(JSON.stringify(suites).length, name).toBeLessThan(1_000_000);
  }
  expect(regressions).toBe(24);
  expect(() => nativeInputs('missing')).toThrow('Unregistered');
});

test('native random-stream and polynomial inputs change with the replay seed', () => {
  for (const suffix of ['ZZ_random_stream.native.json.gz', 'ZZX1_division.native.json.gz']) {
    const name = names.find(n => n.endsWith('/' + suffix))!;
    expect(payloads(name, 123)).not.toEqual(payloads(name, 456));
  }
});
