import { nativeFixtures as loadLiveNative } from "../../../../../tests/property/native-live.mjs";


import { expect, test } from 'bun:test';
import { execFileSync } from 'node:child_process';
const fixtures = (await loadLiveNative(import.meta.url, "./conjugate_bound.native.json.gz")) as { args: (string)[]; error: null; errorType: null; function: string; result: string; seed: number }[];
import { numberofconjugates } from '@sagemath-ts/parigp-ts/src/galconj.js';
import { Flx_deriv, Flx_is_squarefree } from '@sagemath-ts/parigp-ts/src/Flx.js';
import { FpX_red } from '@sagemath-ts/parigp-ts/src/ffinit.js';
import { getrand, setrand } from '@sagemath-ts/parigp-ts/src/random.js';
import { numberofconjugates as legacyCount } from './pari_nf.js';

for (const row of fixtures) {
  test(`native conjugate bound ${row.seed}`, () => {
    const [op, p, start] = row.args.slice(0, 3).map(BigInt);
    const text = row.args[3]!.slice(1, -1).trim();
    const coefficients = text ? text.split(',').map((c) => BigInt(c.trim())) : [];
    const before = [...coefficients],
      saved = getrand();
    try {
      setrand(1n);
      let result: string;
      if (op === 1n && start! < -2n) {
        // A regression to the old signed scan would need 2^63 increments.
        // Timeout is a failure, never an expected native comparison result.
        const program = `import { numberofconjugates } from ${JSON.stringify(new URL('./pari_nf.ts', import.meta.url).pathname)};
          import {getrand,setrand} from '@sagemath-ts/parigp-ts/src/random.js';
          setrand(1n);const value=numberofconjugates(${JSON.stringify(coefficients.map(String))}.map(BigInt),BigInt(${JSON.stringify(String(start))}));
          console.log(JSON.stringify([value,getrand()],(_,v)=>typeof v==='bigint'?String(v):v));`;
        result = execFileSync(process.execPath, ['-e', program], {
          encoding: 'utf8',
          timeout: 10000,
        }).trim();
      } else {
        const value =
          op === 0n
            ? numberofconjugates(coefficients, start!)
            : op === 1n
              ? legacyCount(coefficients, start!)
              : op === 2n
                ? Flx_deriv(FpX_red(coefficients, p!), p!)
                : Flx_is_squarefree(FpX_red(coefficients, p!), p!)
                  ? 1
                  : 0;
        result = JSON.stringify([value, getrand()], (_, v) =>
          typeof v === 'bigint' || typeof v === 'number' ? String(v) : v
        );
      }
      expect({ result, error: null, errorType: null }).toEqual({
        result: row.result,
        error: row.error,
        errorType: row.errorType,
      });
      expect(coefficients).toEqual(before);
    } finally {
      setrand(saved);
    }
  }, 30000);
}
