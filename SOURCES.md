# Source References

This project ports algorithms from the following upstream repositories. The `reference/` directory (gitignored) contains clones of these repos for development reference.

## Repository Versions

These are the commits used as reference for the port. Run `./scripts/clone-references.sh` to populate `reference/`.

| Repository | Commit | Date | URL |
|------------|--------|------|-----|
| sagemath/sage | `28a7d042` | 2026-01-25 | https://github.com/sagemath/sage |
| sagemath/cypari2 | `0a3e42e8` | 2025-12-12 | https://github.com/sagemath/cypari2 |
| pari/pari | `d14c2a9e` | 2026-01-28 | https://pari.math.u-bordeaux.fr/git/pari.git |
| flintlib/flint | `27bdabed` | 2026-01-24 | https://github.com/flintlib/flint |
| libntl/ntl | `be43be35` | 2025-11-07 | https://github.com/libntl/ntl |

## Key Source Paths

| sagemath-ts Module | Primary Source |
|--------------------|----------------|
| `arith/` | `sage/src/sage/arith/misc.py` |
| `rings/integer/` | `sage/src/sage/rings/integer.pyx`, `integer_ring.py` |
| `rings/finite_rings/` | `sage/src/sage/rings/finite_rings/` |
| `rings/polynomial/` | `sage/src/sage/rings/polynomial/` |
| `schemes/elliptic_curves/` | `sage/src/sage/schemes/elliptic_curves/` |
| `parigp-ts` | `pari/src/basemath/`, `cypari2/cypari2/` |

## Updating References

To update to newer commits:

```bash
cd reference/<repo>
git fetch origin
git checkout <commit-or-tag>
```

Then update this file with the new commit hash.

## Version Notes

- **SageMath 10.9.beta4**: Current reference version
- Shallow clones (`--depth 1`) are used by default to save space
- For full history, re-clone without `--depth 1`
