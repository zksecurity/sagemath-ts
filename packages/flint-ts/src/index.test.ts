/**
 * @module flint-ts/index.test
 * @description Smoke test for the package barrel.
 *
 * Regression guard: `index.ts` used to re-export the *interfaces* `nmod_t`,
 * `fmpz_factor` and `fmpz_poly_factor` through a value `export { ... }` clause.
 * TypeScript elides those at compile time, but Bun/Node resolve the barrel at
 * runtime and threw `export 'nmod_t' not found in './nmod_poly.js'`, making
 * `import ... from '@sagemath-ts/flint-ts'` fail outright. `tsc --noEmit`
 * does not catch this, so it needs a runtime test.
 */

import { describe, expect, it } from 'bun:test';

describe('@sagemath-ts/flint-ts barrel', () => {
  it('can be imported at runtime', async () => {
    const mod = await import('./index.js');
    expect(typeof mod).toBe('object');
  });

  it('exposes the documented value exports', async () => {
    const mod = (await import('./index.js')) as Record<string, unknown>;
    const expected = [
      'fmpz',
      'fmpz_init',
      'fmpz_gcd',
      'fmpz_poly',
      'fmpz_poly_init',
      'nmod_poly',
      'nmod_poly_init',
      'nmod_poly_gcd',
    ];
    const missing = expected.filter((n) => mod[n] === undefined);
    expect(missing).toEqual([]);
  });

  it('does not re-export interfaces as runtime values', async () => {
    const mod = (await import('./index.js')) as Record<string, unknown>;
    // These are `interface`s in their defining modules; they must be absent at
    // runtime rather than present-but-undefined or throwing on import.
    for (const n of ['nmod_t', 'fmpz_factor', 'fmpz_poly_factor']) {
      expect(Object.hasOwn(mod, n)).toBe(false);
    }
  });
});
