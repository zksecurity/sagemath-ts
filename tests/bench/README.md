# Benchmarks

This directory provides a simple, reproducible benchmark harness to compare
SageMath vs sagemath-ts. Inputs are generated once and reused by both runtimes
to avoid RNG or generator mismatches.

## Layout

```
tests/bench/
|-- cases/            # Benchmark case definitions (*.bench.json)
|-- inputs/           # Generated inputs (*.bench.inputs.json)
|-- python/           # SageMath runner
|-- typescript/       # sagemath-ts runner (Bun)
`-- results/          # Generated results (gitignored)
```

## Workflow

1) Generate deterministic inputs from case definitions:

```bash
bun run bench:generate
```

2) Run SageMath benchmarks:

```bash
bun run bench:sage -- --input tests/bench/inputs/arith.bench.inputs.json \
  > tests/bench/results/sage.arith.json
```

3) Run sagemath-ts benchmarks:

```bash
bun run bench:ts -- --input tests/bench/inputs/arith.bench.inputs.json \
  > tests/bench/results/ts.arith.json
```

## Case Format

`tests/bench/cases/*.bench.json`:

```json
{
  "module": "arith",
  "cases": [
    {
      "id": "gcd-64bit",
      "function": "gcd",
      "seeds": [1],
      "samples": 200,
      "argGenerators": ["randomBigint(1, 1000)", "randomBigint(1, 1000)"],
      "warmup": 1000,
      "iterations": 10000,
      "repeats": 5
    }
  ]
}
```

- `samples`: number of distinct argument tuples generated per seed
- `warmup`: warmup iterations before timing
- `iterations`: timed loop iterations per repeat
- `repeats`: number of timed repeats

Additional generators supported by `tests/bench/generate.ts`:
- `randomBigintBits(bits)` for fixed-width integers (e.g. `randomBigintBits(2048)` for RSA-size inputs).

## Input Format

`tests/bench/inputs/*.bench.inputs.json`:

```json
{
  "module": "arith",
  "generatedAt": "2026-01-30T00:00:00.000Z",
  "sourceFile": "tests/bench/cases/arith.bench.json",
  "cases": [
    {
      "id": "gcd-64bit",
      "function": "gcd",
      "inputs": [["12", "8"], ["99", "33"]],
      "warmup": 1000,
      "iterations": 10000,
      "repeats": 5
    }
  ]
}
```

## Results

Both runners output JSON with total and per-call timing stats. Compare results
by looking at `per_call_ns` (median and p90) across the two outputs.

Notes:
- The benchmark suite includes RSA-2048-style modular exponentiation and
  256-bit curve scalar multiplication (secp256k1, P-256, BLS12-381 G1).
- Pairing benchmarks for BLS12-381 are not included yet because sagemath-ts
  currently only constructs prime fields via `GF(p)` (no extension fields).
