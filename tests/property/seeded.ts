import { createHash, randomBytes } from 'node:crypto';
import type { CaseSuite, CaseRow, CaseValue } from './case-format.js';
import { decodeArg, encodeArg, isRowCase, parseFixedGenerator } from './case-format.js';
import { MersenneTwister } from './typescript/mersenne-twister.js';

/** Version the input algorithm independently of the implementation under test. */
export const GENERATOR_VERSION = 'mt19937-recipes-v2';

export function freshSeed(): number {
  return randomBytes(4).readUInt32LE();
}

export function derivedSeed(seed: number, label: string): number {
  return createHash('sha256').update(`${GENERATOR_VERSION}:${seed}:${label}`).digest().readUInt32LE();
}

// A test-input helper, deliberately independent of the port's primality routine.
function prime(n: bigint): boolean {
  if (n < 2n) return false;
  for (let d = 2n; d * d <= n; d++) if (n % d === 0n) return false;
  return true;
}

export function generateArgument(spec: string, rng: MersenneTwister): CaseValue {
  const fixed = parseFixedGenerator(spec);
  if (fixed !== undefined) return fixed;
  const integer = /^random(Bigint|Prime)\((-?\d+),\s*(-?\d+)\)$/.exec(spec);
  if (integer) {
    const lo = BigInt(integer[2]!), hi = BigInt(integer[3]!);
    if (lo > hi) throw new Error(`Empty generator interval: ${spec}`);
    if (integer[1] === 'Bigint') return rng.randint(lo, hi);
    if (hi > 1_000_000n || lo < 0n) throw new Error(`Prime generator requires 0 <= min <= max <= 1000000: ${spec}`);
    const values: bigint[] = [];
    for (let n = lo; n <= hi; n++) if (prime(n)) values.push(n);
    if (!values.length) throw new Error(`No prime in generator interval: ${spec}`);
    return values[Number(rng.randint(0n, BigInt(values.length - 1)))]!;
  }
  if (spec.startsWith('randomList(') && spec.endsWith(')')) {
    const body = spec.slice(11, -1), comma = body.lastIndexOf(',');
    const length = Number(body.slice(comma + 1).trim());
    if (!Number.isSafeInteger(length) || length < 0 || length > 100_000)
      throw new Error(`Invalid list length: ${spec}`);
    return Array.from({ length }, () => generateArgument(body.slice(0, comma).trim(), rng)).flat();
  }
  throw new Error(`Unknown argument generator: ${spec}`);
}

export type ArgumentRecipe =
  | { constant: import('./case-format.js').JsonArg }
  | { integer: [number | string, number | string] }
  | { bytes: number }
  | { shifted: { coefficient: number; shift: number; offset: number }[] }
  | { text: 'expression' | 'identifier' }
  | { list: { length: number; range: [number, number]; nonzeroLast?: boolean } }
  | { scale: { values: (number | string)[]; range: [number, number] } };

export function generateRecipe(recipe: ArgumentRecipe[], rng: MersenneTwister): CaseValue[] {
  const integer = ([lo, hi]: [number | string, number | string]) => {
    if (BigInt(lo) > BigInt(hi)) throw new Error('Empty recipe interval');
    return rng.randint(BigInt(lo), BigInt(hi));
  };
  const list = (length: number, value: () => bigint) => {
    if (!Number.isSafeInteger(length) || length < 0 || length > 100_000)
      throw new Error('Invalid recipe list length');
    return Array.from({ length }, value);
  };
  return recipe.map(spec => {
    if ('shifted' in spec) return spec.shifted.map(v => {
      if (!Number.isSafeInteger(v.shift) || v.shift < 0 || v.shift > 100_000) throw new Error('Invalid shift');
      return (BigInt(v.coefficient) << BigInt(v.shift)) + BigInt(v.offset);
    });
    if ('text' in spec) {
      const a = rng.randint(0n, 10000n), b = rng.randint(1n, 10000n);
      const text = spec.text === 'expression' ? `(${a} + ${b}) * ${rng.randint(0n, 100n)}` : `v${a}_${b}`;
      return Array.from(text, c => BigInt(c.charCodeAt(0)));
    }
    if ('constant' in spec) return decodeArg(spec.constant);
    if ('integer' in spec) return integer(spec.integer);
    if ('bytes' in spec) return list(spec.bytes, () => rng.randint(0n, 255n));
    if ('list' in spec) {
      const values = list(spec.list.length, () => integer(spec.list.range));
      if (spec.list.nonzeroLast) {
        const [lo, hi] = spec.list.range;
        if (!values.length || (lo === 0 && hi === 0)) throw new Error('Nonzero final coefficient requires a nonempty domain');
        // Uniformly omit zero without an unbounded rejection loop.
        const hasZero = lo <= 0 && hi >= 0;
        let last = integer([lo, hi - (hasZero ? 1 : 0)]);
        if (hasZero && last >= 0n) last++;
        values[values.length - 1] = last;
      }
      return values;
    }
    if ('scale' in spec) {
      const factor = integer(spec.scale.range);
      return spec.scale.values.map(v => BigInt(v) * factor);
    }
    throw new Error(`Unknown argument recipe: ${JSON.stringify(spec)}`);
  });
}

/** Materialize once so both runtimes receive byte-identical, lossless inputs. */
export function materializeSuite(suite: CaseSuite, seed: number, runs: number): CaseSuite {
  if (!Number.isSafeInteger(seed) || seed < 0 || seed > 0xffffffff)
    throw new Error('seed must be an unsigned 32-bit integer');
  if (!Number.isSafeInteger(runs) || runs < 1 || runs > 100_000)
    throw new Error('runs must be an integer between 1 and 100000');
  return {
    module: suite.module,
    cases: suite.cases.map((c, index) => {
      if (isRowCase(c)) return c;
      const specs = c.argGenerators ?? [];
      const recipes = c.recipes as ArgumentRecipe[][] | undefined;
      if (recipes && !recipes.length) throw new Error(`No recipes for ${c.function}`);
      const rows: CaseRow[] = [];
      const replays = (c.replaySeeds ?? []) as { seed: number; runs: number; reason: string }[];
      const seen = new Set<string>();
      for (const [planIndex, plan] of [{ seed, runs }, ...replays].entries()) {
        materializeSuite({ module: 'validation', cases: [] }, plan.seed, plan.runs);
        const seeds = [
          ...(planIndex === 0 ? c.seeds ?? [] : []),
          ...Array.from({ length: Math.max(plan.runs, recipes?.length ?? 0) }, (_, trial) =>
            derivedSeed(plan.seed, `${suite.module}:${c.function}:${index}:${trial}`)),
        ];
        for (const [trial, caseSeed] of seeds.entries()) {
          const rng = new MersenneTwister(caseSeed);
          // Visit each constructor variant and replay retained failures alongside
          // fresh trials. Recipe/index changes require reviewing retained seeds.
          const values = recipes
            ? generateRecipe(recipes[(trial + derivedSeed(plan.seed, `${suite.module}:${c.function}:variant`)) % recipes.length]!, rng)
            : specs.map(spec => generateArgument(spec, rng));
          const row: CaseRow = [caseSeed, ...values.map(encodeArg)];
          const key = JSON.stringify(row);
          if (!seen.has(key)) { seen.add(key); rows.push(row); }
        }
      }
      return { function: c.function, rows };
    }),
  };
}
