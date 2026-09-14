/** 64-bit word random kernels from ulong_extras.h and ulong_extras/randomisation.c. */
import type { flint_rand_t } from '../flint.js';
export function n_randlimb(state: flint_rand_t): bigint {
  state.__randval = BigInt.asUintN(64, state.__randval * 13282407956253574709n + 286824421n);
  state.__randval2 = BigInt.asUintN(64, state.__randval2 * 7557322358563246341n + 286824421n);
  return (state.__randval >> 32n) + ((state.__randval2 >> 32n) << 32n);
}
export function n_randint(state: flint_rand_t, limit: bigint): bigint {
  const n = BigInt.asUintN(64, limit),
    r = n_randlimb(state);
  return (n & (n - 1n)) === 0n ? r & (n - 1n) : (r * n) >> 64n;
}
export function n_urandint(state: flint_rand_t, limit: bigint): bigint {
  return n_randint(state, limit);
}
/** @see Deviation: FLINT Random Bit-Count Bounds */
export function n_randbits(state: flint_rand_t, bits: number): bigint {
  if (!Number.isInteger(bits) || bits < 0 || bits > 64)
    throw new RangeError('bits must be between 0 and 64');
  if (bits === 0) return 0n;
  return (1n << BigInt(bits - 1)) | n_randint(state, bits === 64 ? 0n : 1n << BigInt(bits));
}
/** @see Deviation: FLINT Random Bit-Count Bounds */
export function n_randtest_bits(state: flint_rand_t, bits: number): bigint {
  if (!Number.isInteger(bits) || bits < 0 || bits > 64)
    throw new RangeError('bits must be between 0 and 64');
  const m = n_randlimb(state);
  if ((m & 7n) !== 0n) return n_randbits(state, bits);
  let n: bigint;
  switch (Number((m >> 3n) & 7n)) {
    case 0:
      n = 0n;
      break;
    case 1:
      n = 1n;
      break;
    case 2:
      n = (1n << 62n) - 1n;
      break;
    case 3:
      n = (1n << 63n) - 1n;
      break;
    case 4:
      n = (1n << 64n) - 1n;
      break;
    case 5:
      n = (1n << n_randint(state, 64n)) - (1n << n_randint(state, 64n));
      break;
    case 6:
      n = 1n << n_randint(state, 64n);
      break;
    default:
      n = -(1n << n_randint(state, 64n));
  }
  n = BigInt.asUintN(64, n);
  if (bits < 64) n &= (1n << BigInt(bits)) - 1n;
  return bits ? n | (1n << BigInt(bits - 1)) : 0n;
}
export function n_randtest(state: flint_rand_t): bigint {
  return n_randtest_bits(state, Number(n_randint(state, 65n)));
}
export function n_randtest_not_zero(state: flint_rand_t): bigint {
  let n: bigint;
  do {
    n = n_randtest(state);
  } while (n === 0n);
  return n;
}
