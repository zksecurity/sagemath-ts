# sagemath-ts

> ⚠️ **EXPERIMENTAL - DO NOT USE IN PRODUCTION**
>
> This library is **experimental** and provided for educational and research purposes only.
>
> - **Not production-ready**: Code has not been audited and may contain bugs
> - **Not constant-time**: Implementations are vulnerable to timing side-channel attacks
> - **AI-assisted development**: Significant portions were written by AI agents
> - **Partial port**: Only number theory modules are implemented (not a complete SageMath replacement)
> - **No security guarantees**: Do not use for real cryptographic applications
>
> For production cryptography, use established libraries like [OpenSSL](https://www.openssl.org/), [libsodium](https://libsodium.org/), or [noble-curves](https://github.com/paulmillr/noble-curves).

A TypeScript port of the **number theory modules** from [SageMath](https://github.com/sagemath/sage), focusing on educational exploration of cryptographic primitives. It's ~3.5MB, as opposed to the 5GB of SageMath when installed as a Linux package!

## Scope

This project ports **only the number theory and algebraic modules** from SageMath—specifically those relevant to cryptography education. It is **not** a complete SageMath replacement. Modules like symbolic computation, calculus, plotting, and statistics are not included.

## Goals

1. **Mirror SageMath structure** - Same file paths, function names, and module organization
2. **Deterministic equivalence** - Property tests ensure identical outputs for same inputs
3. **Cryptography education** - Prioritize number theory modules used in cryptographic applications

## ⚠️ Behavioral Differences

This port aims for **exact behavioral equivalence** with SageMath, but some differences are intentional. **Before using this library, review [`DEVIATIONS.md`](./DEVIATIONS.md)** which documents:

- Algorithm implementation choices (e.g., Stein's GCD vs Euclidean)
- Error handling differences (`null` returns vs exceptions)
- Type system adaptations for TypeScript
- Feature limitations (e.g., tower field extension degrees)

Every deviation includes rationale and behavioral impact.

## Project Structure

```
sagemath-ts/
├── reference/              # Cloned upstream repos for reference
│   ├── sage/              # SageMath main repository
│   ├── cypari2/           # Python PARI/GP interface
│   ├── pari/              # PARI/GP source
│   ├── flint/             # FLINT library
│   └── ntl/               # NTL library
├── packages/              # TypeScript implementations
│   ├── sagemath-ts/       # Main package (mirrors sage/src/sage/)
│   ├── parigp-ts/         # PARI/GP bindings port
│   ├── flint-ts/          # FLINT bindings port
│   ├── ntl-ts/            # NTL bindings port
│   └── zksecurity-cheatsheets/ # Crypto parameter cheatsheets
├── tests/
│   ├── property/          # Cross-language property tests
│   └── unit/              # TypeScript unit tests
└── docs/                  # Documentation and style guides
```

## Setup

```bash
# Clone reference repositories
./scripts/clone-references.sh

# Install dependencies
bun install

# Run tests
bun test
```

## Playground

An interactive browser playground lets you experiment with sagemath-ts in real time, with 77 tutorial lessons covering number theory and cryptography.

```bash
# Start the playground server (builds automatically if needed)
bun playground

# Start with hot reload (watches for file changes and rebuilds)
bun playground:dev
```

Open http://localhost:5173 in your browser. The server serves:
- **Lessons** - 77 interactive tutorials from `tutorial/` with runnable code
- **API Docs** - Generated documentation for all exported functions
- **Playground** - Free-form code editor for experimentation

To rebuild manually:

```bash
bun playground/build.ts
```

Lessons are auto-generated from the TypeScript files in `tutorial/`. Edit those files and the playground will update (with hot reload enabled).

## Property Testing Strategy

We ensure deterministic equivalence between Python/SageMath and TypeScript implementations:

1. **Seed both environments** with identical random seeds
2. **Execute same operations** via parallel test files (`.py` and `.ts`)
3. **Compare transcripts** - outputs must match exactly

See `tests/property/README.md` for details.

## Benchmarks (SageMath vs sagemath-ts)

These are microbenchmarks using identical, pre-generated inputs. Results are from
`tests/bench/results/*.json` generated on **January 30, 2026** using **Bun 1.3.6**
and **Sage 10.3**.

Median per-call timings (p50, microseconds) and relative speed (TS / Sage):

| Case | TS p50 (µs) | Sage p50 (µs) | TS / Sage |
| --- | ---: | ---: | ---: |
| gcd-64bit | 8.14 | 1.18 | 6.9× slower |
| is-prime-32bit | 1.57 | 1.18 | 1.3× slower |
| factor-32bit | 66.31 | 5.67 | 11.7× slower |
| power-mod-32bit | 1.58 | 4.98 | 0.32× (≈3.2× faster) |

Crypto-scale microbenchmarks (same machine/date; inputs in `tests/bench/inputs/crypto.bench.inputs.json`):

| Case | TS p50 (µs) | Sage p50 (µs) | TS / Sage |
| --- | ---: | ---: | ---: |
| rsa-2048-pow | 54.81 | 24.94 | 2.20× slower |
| secp256k1-mul | 799.59 | 322.97 | 2.48× slower |
| p256-mul | 841.91 | 345.12 | 2.44× slower |
| bls12-381-g1-mul | 1129.17 | 476.51 | 2.37× slower |

Notes:
- These are **not** end-to-end benchmarks; they measure tight loops over small inputs.
- The RSA-2048 case is a pure modular exponentiation with random 2048-bit moduli (not full keygen/sign/verify).
- EC scalar multiplication uses the generic short-Weierstrass formulas in this repo (not specialized, side-channel-hardened code).
- Expect different numbers across machines and runtimes; rerun with `bun run bench:generate`,
  then `bun run bench:ts` / `bun run bench:sage` if you need local data.

## Documentation

| Document                                                 | Purpose                                |
| -------------------------------------------------------- | -------------------------------------- |
| [`SCOPE.md`](./SCOPE.md)                                 | Implementation progress tracker        |
| [`SOURCES.md`](./SOURCES.md)                             | Upstream repo commits used for porting |
| [`DEVIATIONS.md`](./DEVIATIONS.md)                       | Behavioral differences from SageMath   |
| [`AGENTS.md`](./AGENTS.md)                               | Guidelines for AI agents               |
| [`docs/PORTING_GUIDE.md`](./docs/PORTING_GUIDE.md)       | Standards for porting new modules      |
| [`tests/property/README.md`](./tests/property/README.md) | Property testing framework             |

## Contributing

1. Read the [Porting Guide](./docs/PORTING_GUIDE.md) for standards
2. Check [SCOPE.md](./SCOPE.md) for available work
3. Follow [AGENTS.md](./AGENTS.md) guidelines
4. Document any deviations in [DEVIATIONS.md](./DEVIATIONS.md)

## License

Same license as SageMath (GPL-2.0+).
