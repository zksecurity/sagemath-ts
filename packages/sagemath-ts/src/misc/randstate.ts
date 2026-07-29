/**
 * @module sage/misc/randstate
 * @description Random state utilities
 *
 * Port of: sage/misc/randstate.pyx
 * Reference: reference/sage/src/sage/misc/randstate.pyx
 *
 * Sage's `randstate` wraps a GMP `gmp_randstate_t` initialised with
 * `gmp_randinit_default` (`randstate.pyx:511`), which is GMP's Mersenne
 * Twister (`gmp_randinit_mt`), and seeded with `gmp_randseed`
 * (`randstate.pyx:558`).  All raw bits are drawn through
 * `gmp_urandomb_ui(state, nbits)` (`randstate.pyx:881`).
 *
 * Everything here is a port of the real upstream sources, so a seeded stream
 * reproduces Sage's bit-for-bit:
 *
 * - the generator and its buffer bookkeeping: GMP 6.3.0
 *   `rand/randmt.c` (`__gmp_mt_recalc_buffer`, `__gmp_randget_mt`,
 *   `__gmp_randinit_mt_noseed`);
 * - the seeding permutation: GMP 6.3.0 `rand/randmts.c`
 *   (`mangle_seed`, `randseed_mt`);
 * - `random_below`: GMP 6.3.0 `mpz/urandomm.c` (`mpz_urandomm`);
 * - {@link PythonRandom}: CPython's `Modules/_randommodule.c` plus
 *   `Lib/random.py`, which is the *second* generator Sage exposes
 *   (`randstate.pyx:583-627` `randstate.python_random`) and the one that
 *   `sage.misc.prandom` — hence `IntegerModRing.random_element`
 *   (`integer_mod_ring.py:1543-1547`) — actually draws from.
 */

import { ValueError } from '../errors.js';

// Sage uses a 31-bit RAND_MAX for some distributions (`randstate.pyx:60`).
export const SAGE_RAND_MAX = 0x7fffffff;

/** Number of 32-bit words in the MT19937 state (`randmt.h:36`). */
const MT_N = 624;
/** MT19937 shift parameter (`randmt.h:37`). */
const MT_M = 397;
/** MT19937 twist constant (`randmt.h:38`). */
const MT_MATRIX_A = 0x9908b0df;
/** Most significant bit. */
const MT_UPPER_MASK = 0x80000000;
/** Least significant 31 bits. */
const MT_LOWER_MASK = 0x7fffffff;
/** Number of extractions used to warm the buffer up (`randmt.h:33`). */
const MT_WARM_UP = 2000;

/**
 * Iteration limit in `mpz_urandomm` (`gmp-impl.h`, `MAX_URANDOMM_ITER`).
 */
const MAX_URANDOMM_ITER = 80;

/**
 * GMP's initial Mersenne Twister buffer for an *unseeded* generator.
 *
 * Verbatim from `default_state[N]` in GMP 6.3.0 `rand/randmt.c:57-176`.  Sage
 * reaches this state whenever the seed is `0`, because `randstate.__init__`
 * skips `gmp_randseed` in that case (`randstate.pyx:556-559`).
 */
const MT_DEFAULT_STATE: readonly number[] = [
  0xd247b233, 0x9e5aa8f1, 0x0ffa981b, 0x9dcb0980, 0x74200f2b, 0xa576d044, 0xe9f05adf, 0x1538bff5,
  0x59818bbf, 0xcf9e58d8, 0x09fce032, 0x6a1c663f, 0x5116e78a, 0x69b3e0fa, 0x6d92d665, 0xd0a8be98,
  0xf669b734, 0x41ac1b68, 0x630423f1, 0x4b8d6b8a, 0xc2c46dd7, 0x5680747d, 0x43703e8f, 0x3b6103d2,
  0x49e5eb3f, 0xcbdab4c1, 0x9c988e23, 0x747bee0b, 0x9111e329, 0x9f031b5a, 0xecca71b9, 0x2afe4ef8,
  0x8421c7ed, 0xac89aff1, 0xaed90df3, 0x2dd74f01, 0x14906a13, 0x75873fa9, 0xff83f877, 0x5028a0c9,
  0x11b4c41d, 0x7caedbc4, 0x8672d0a7, 0x48a7c109, 0x8320e59f, 0xbc0b3d5f, 0x75a30886, 0xf9e0d128,
  0x41af7580, 0x239bb94d, 0xc67a3c81, 0x74eebd6e, 0xbc02b53c, 0x727ea449, 0x6b8a2806, 0x5853b0da,
  0xbde032f4, 0xce234885, 0x320d6145, 0x48cc053f, 0x00dbc4d2, 0xd55a2397, 0xe1059b6f, 0x1c3e05d1,
  0x09657c64, 0xd07cb661, 0x6e982e34, 0x6dd1d777, 0xeded1071, 0xd79dfd65, 0xf816ddce, 0xb6faf1e4,
  0x1c771074, 0x311835bd, 0x18f952f7, 0xf8f40350, 0x4eced354, 0x7c8ac12b, 0x31a9994d, 0x4fd47747,
  0xdc227a23, 0x6dfafddf, 0x6796e748, 0x0c6f634f, 0xf992fa1d, 0x4cf670c9, 0x067dfd31, 0xa7a3e1a5,
  0x8cd7d9df, 0x972ccb34, 0x67c82156, 0xd548f6a8, 0x045cec21, 0xf3240bfb, 0xdef656a7, 0x43de08c5,
  0xdad1f92f, 0x3726c56b, 0x1409f19a, 0x942fd147, 0xb926749c, 0xaddc31b8, 0x53d0d869, 0xd1ba52fe,
  0x6722df8c, 0x22d95a74, 0x7dc1b52a, 0x1dec6fd5, 0x7262874d, 0x0a725dc9, 0xe6a8193d, 0xa052835a,
  0xdc9ad928, 0xe59ebb90, 0x70dba9ff, 0xd612749d, 0x5a5a638c, 0x6086ec37, 0x2a579709, 0x1449ea3a,
  0xbc8e3c06, 0x2f900666, 0xfbe74fd1, 0x6b35b911, 0xf8335008, 0xef1e979d, 0x738ab29d, 0xa2dc0fdc,
  0x7696305d, 0xf5429dac, 0x8c41813b, 0x8073e02e, 0xbef83ccd, 0x7b50a95a, 0x05ee5862, 0x00829ece,
  0x8ca1958c, 0xbe4ea2e2, 0x4293bb73, 0x656f7b23, 0x417316d8, 0x4467d7cf, 0x2200e63b, 0x109050c8,
  0x814cbe47, 0x36b1d4a8, 0x36af9305, 0x308327b3, 0xebcd7344, 0xa738de27, 0x5a10c399, 0x4142371d,
  0x64a18528, 0x0b31e8b2, 0x641057b9, 0x6afc363b, 0x108ad953, 0x9d4da234, 0x0c2d9159, 0x1c8a1a1f,
  0x310c66ba, 0x87aa1070, 0xdac832ff, 0x0a433422, 0x7af15812, 0x2d8d9bd0, 0x995a25e9, 0x25326cac,
  0xa34384db, 0x4c8421cc, 0x4f0315ec, 0x29e8649e, 0xa7732d6f, 0x2e94d3e3, 0x7d98a340, 0x397c4d74,
  0x659db4de, 0x747d4e9a, 0xd9db8435, 0x4659dbe9, 0x313e6dc5, 0x29d104dc, 0x9f226cba, 0x452f18b0,
  0xd0bc5068, 0x844ca299, 0x782b294e, 0x4ae2eb7b, 0xa4c475f8, 0x70a81311, 0x4b3e8bcc, 0x7e20d4ba,
  0xabca33c9, 0x57be2960, 0x44f9b419, 0x2e567746, 0x72eb757a, 0x102cc0e8, 0xb07f32b9, 0xd0dabd59,
  0xba85ad6b, 0xf3e20667, 0x98d77d81, 0x197afa47, 0x518ee9ac, 0xe10ce5a2, 0x01cf2c2a, 0xd3a3af3d,
  0x16ddfd65, 0x669232f8, 0x1c50a301, 0xb93d9151, 0x9354d3f4, 0x847d79d0, 0xd5fe2ec6, 0x1f7b0610,
  0xfa6b90a5, 0xc5879041, 0x2e7dc05e, 0x423f1f32, 0xef623ddb, 0x49c13280, 0x98714e92, 0xc7b6e4ad,
  0xc4318466, 0x0737f312, 0x4d3c003f, 0x9acc1f1f, 0x5f1c926d, 0x085fa771, 0x185a83a2, 0xf9aa159d,
  0x0b0b0132, 0xf98e7a43, 0xcd9ebdbe, 0x0190cb29, 0x10d93fb6, 0x3b8a4d97, 0x66a65a41, 0xe43e766f,
  0x77be3c41, 0xb9686364, 0xcb36994d, 0x6846a287, 0x567e77f7, 0x36178dd8, 0xbde6b1f2, 0xb6efdc64,
  0x82950324, 0x42053f47, 0xc09be51c, 0x0942d762, 0x35f92c7f, 0x367dec61, 0x6ee3d983, 0xdbaaf78a,
  0x265d2c47, 0x8eb4bf5c, 0x33b232d7, 0xb0137e77, 0x373c39a7, 0x8d2b2e76, 0xc7510f01, 0x50f9e032,
  0x7b1fdddb, 0x724c2aae, 0xb10ecb31, 0xcca3d1b8, 0x7f0bcf10, 0x4254bbbd, 0xe3f93b97, 0x2305039b,
  0x53120e22, 0x1a2f3b9a, 0x0fddbd97, 0x0118561e, 0x0a798e13, 0x9e0b3acd, 0xdb6c9f15, 0xf512d0a2,
  0x9e8c3a28, 0xee2184ae, 0x0051ec2f, 0x2432f74f, 0xb0aa66ea, 0x55128d88, 0xf7d83a38, 0x4dae8e82,
  0x3fdc98d6, 0x5f0bd341, 0x7244be1d, 0xc7b48e78, 0x2d473053, 0x43892e20, 0xba0f1f2a, 0x524d4895,
  0x2e10bcb1, 0x4c372d81, 0x5c3e50cd, 0xcf61cc2e, 0x931709ab, 0x81b3aefc, 0x39e9405e, 0x7ffe108c,
  0x4fbb3ff8, 0x06abe450, 0x7f5bf51e, 0xa4e3cdfd, 0xdb0f6c6f, 0x159a1227, 0x3b9fed55, 0xd20b6f7f,
  0xfbe9cc83, 0x64856619, 0xbf52b8af, 0x9d7006b0, 0x71165bc6, 0xae324aee, 0x29d27f2c, 0x794c2086,
  0x74445ce2, 0x782915cc, 0xd4ce6886, 0x3289ae7c, 0x53def297, 0x4185f7ed, 0x88b72400, 0x3c09dc11,
  0xbce3aab6, 0x6a75934a, 0xb267e399, 0x000df1bf, 0x193ba5e2, 0xfa3e1977, 0x179e14f6, 0x1eede298,
  0x691f0b06, 0xb84f78ac, 0xc1c15316, 0xffff3ad6, 0x0b457383, 0x518cd612, 0x05a00f3e, 0xd5b7d275,
  0x4c5eccd7, 0xe02cd0be, 0x5558e9f2, 0x0c89bbf0, 0xa3d96227, 0x2832d2b2, 0xf667b897, 0xd4556554,
  0xf9d2f01f, 0xfa1e3fae, 0x52c2e1ee, 0xe5451f31, 0x7e849729, 0xdabdb67a, 0x54bf5e7e, 0xf831c271,
  0x5f1a17e3, 0x9d140afe, 0x92741c47, 0x48cfabce, 0x9cbbe477, 0x9c3ee57f, 0xb07d4c39, 0xcc21bce2,
  0x697708b1, 0x58da2a6b, 0x2370db16, 0x6e641948, 0xacc5bd52, 0x868f24cc, 0xca1db0f5, 0x4cada492,
  0x3f443e54, 0xc4a4d5e9, 0xf00ad670, 0xe93c86e0, 0xfe90651a, 0xdde532a3, 0xa66458df, 0xab7d7151,
  0x0e2e775f, 0xc9109f99, 0x8d96d59f, 0x73cef14c, 0xc74e88e9, 0x02712dc0, 0x04f41735, 0x2e5914a2,
  0x59f4b2fb, 0x0287fc83, 0x80bc0343, 0xf6b32559, 0xc74178d4, 0xf1d99123, 0x383ccc07, 0xacc0637d,
  0x0863a548, 0xa6fcac85, 0x2a13eff0, 0xaf2eedb1, 0x41e72750, 0xe0c6b342, 0x5da22b46, 0x635559e0,
  0xd2ea40ac, 0x10aa98c0, 0x19096497, 0x112c542b, 0x2c85040c, 0xa868e7d0, 0x6e260188, 0xf596d390,
  0xc3bb5d7a, 0x7a2aa937, 0xdfd15032, 0x6780ae3b, 0xdb5f9cd8, 0x8bd266b0, 0x7744af12, 0xb463b1b0,
  0x589629c9, 0xe30dbc6e, 0x880f5569, 0x209e6e16, 0x9deca50c, 0x02987a57, 0xbed3ea57, 0xd3a678aa,
  0x70dd030d, 0x0cfd9c5d, 0x92a18e99, 0xf5740619, 0x7f6f0a7d, 0x134caf9a, 0x70f5bae4, 0x23dca7b5,
  0x4d788fcd, 0xc7f07847, 0xbcf77da1, 0x9071d568, 0xfc627ea1, 0xae004b77, 0x66b54bcb, 0x7ef2daac,
  0xdcd5ac30, 0xb9bdf730, 0x505a97a7, 0x9d881fd3, 0xadb796cc, 0x94a1d202, 0x97535d7f, 0x31ec20c0,
  0xb1887a98, 0xc1475069, 0xa6f73af3, 0x71e4e067, 0x46a569de, 0xd2ade430, 0x6f0762c7, 0xf50876f4,
  0x53510542, 0x03741c3e, 0x53502224, 0xd8e54d60, 0x3c44ab1a, 0x34972b46, 0x74bfa89d, 0xd7d768e0,
  0x37e605dc, 0xe13d1bdf, 0x5051c421, 0xb9e057be, 0xb717a14c, 0xa1730c43, 0xb99638be, 0xb5d5f36d,
  0xe960d9ea, 0x6b1388d3, 0xecb6d3b6, 0xbdbe8b83, 0x2e29afc5, 0x764d71ec, 0x4b8f4f43, 0xc21ddc00,
  0xa63f657f, 0x82678130, 0xdbf535ac, 0xa594fc58, 0x942686bc, 0xbd9b657b, 0x4a0f9b61, 0x44ff184f,
  0x38e10a2f, 0x61910626, 0x5e247636, 0x7106d137, 0xc62802f0, 0xbd1d1f00, 0x7cc0dcb2, 0xed634909,
  0xdc13b24e, 0x9799c499, 0xd77e3d6a, 0x14773b68, 0x967a4fb7, 0x35eecfb1, 0x2a5110b8, 0xe2f0af94,
  0x9d09dea5, 0x20255d27, 0x5771d34b, 0xe1089ee4, 0x246f330b, 0x8f7caee5, 0xd3064712, 0x75cafbee,
  0xb94f7028, 0xed953666, 0x5d1975b4, 0x5af81271, 0x13be2025, 0x85194659, 0x30805331, 0xec9d46c0,
  0xbc027c36, 0x2af84188, 0xc2141b80, 0xc02b1e4a, 0x04d36177, 0xfc50e9d7, 0x39ce79da, 0x917e0a00,
  0xef7a0bf4, 0xa98bd8d1, 0x19424dd2, 0x9439df1f, 0xc42af746, 0xaddbe83e, 0x85221f0d, 0x45563e90,
  0x9095ec52, 0x77887b25, 0x8ae46064, 0xbd43b71a, 0xbb541956, 0x7366cf9d, 0xee8e1737, 0xb5a727c9,
  0x5076b3e7, 0xfc70baca, 0xce135b75, 0xc4e91aa3, 0xf0341911, 0x53430c3f, 0x886b0824, 0x6bb5b8b7,
  0x33e21254, 0xf193b456, 0x5b09617f, 0x215fff50, 0x48d97ef1, 0x356479ab, 0x6ea9ddc4, 0x0d352746,
  0xa2f5ce43, 0xb226a1b3, 0x1329ea3c, 0x7a337cc2, 0xb5cce13d, 0x563e3b5b, 0x534e8e8f, 0x561399c9,
  0xe1596392, 0xb0f03125, 0x4586645b, 0x1f371847, 0x94eaabd1, 0x41f97edd, 0xe3e5a39b, 0x71c774e2,
  0x507296f4, 0x5960133b, 0x7852c494, 0x3f5b2691, 0xa3f87774, 0x5a7af89e, 0x17da3f28, 0xe9d9516d,
  0xfcc1c1d5, 0xe4618628, 0x04081047, 0xd8e4db5f, 0xdc380416, 0x8c4933e2, 0x95074d53, 0xb1b0032d,
  0xcc8102ea, 0x71641243, 0x98d6eb6a, 0x90fec945, 0xa0914345, 0x6fab037d, 0x70f49c4d, 0x05bf5b0e,
  0x927aaf7f, 0xa1940f61, 0xfee0756f, 0xf815369f, 0x5c00253b, 0xf2b9762f, 0x4aeb3ccc, 0x1069f386,
  0xfba4e7b9, 0x70332665, 0x6bca810e, 0x85ab8058, 0xae4b2b2f, 0x9d120712, 0xbee8eacb, 0x776a1112,
];

let _currentRandState: RandState | null = null;

function defaultSeed(): bigint {
  const timeSeed = BigInt(Date.now());
  const perfSeed =
    typeof performance !== 'undefined' ? BigInt(Math.floor(performance.now() * 1_000_000)) : 0n;
  const hrSeed =
    typeof process !== 'undefined' && typeof process.hrtime?.bigint === 'function'
      ? process.hrtime.bigint()
      : 0n;

  if (typeof globalThis.crypto?.getRandomValues === 'function') {
    const buf = new Uint32Array(2);
    globalThis.crypto.getRandomValues(buf);
    const cryptoSeed = (BigInt(buf[0]) << 32n) | BigInt(buf[1]);
    return cryptoSeed ^ timeSeed ^ perfSeed ^ hrSeed;
  }

  return timeSeed ^ perfSeed ^ hrSeed;
}

/**
 * MT19937 buffer regeneration.
 *
 * Verbatim port of `__gmp_mt_recalc_buffer` (GMP `rand/randmt.c:179-198`).
 * CPython's `genrand_uint32` regenerates the buffer with exactly the same
 * loop (`Modules/_randommodule.c`), so both generators share it.
 */
function mtRecalcBuffer(mt: Uint32Array): void {
  let kk = 0;
  let y: number;
  for (; kk < MT_N - MT_M; kk++) {
    y = (mt[kk]! & MT_UPPER_MASK) | (mt[kk + 1]! & MT_LOWER_MASK);
    mt[kk] = (mt[kk + MT_M]! ^ (y >>> 1) ^ (y & 1 ? MT_MATRIX_A : 0)) >>> 0;
  }
  for (; kk < MT_N - 1; kk++) {
    y = (mt[kk]! & MT_UPPER_MASK) | (mt[kk + 1]! & MT_LOWER_MASK);
    mt[kk] = (mt[kk - (MT_N - MT_M)]! ^ (y >>> 1) ^ (y & 1 ? MT_MATRIX_A : 0)) >>> 0;
  }
  y = (mt[MT_N - 1]! & MT_UPPER_MASK) | (mt[0]! & MT_LOWER_MASK);
  mt[MT_N - 1] = (mt[MT_M - 1]! ^ (y >>> 1) ^ (y & 1 ? MT_MATRIX_A : 0)) >>> 0;
}

/**
 * MT19937 tempering, shared by GMP's `NEXT_RANDOM` (`randmt.c:213-224`) and
 * CPython's `genrand_uint32`.
 */
function mtTemper(word: number): number {
  let y = word;
  y ^= y >>> 11;
  y = (y ^ ((y << 7) & 0x9d2c5680)) >>> 0;
  y = (y ^ ((y << 15) & 0xefc60000)) >>> 0;
  y ^= y >>> 18;
  return y >>> 0;
}

/** `2^19937 - 20023`, the powering modulus of `mangle_seed` (`randmts.c:35`). */
const MANGLE_SHIFT = 19937n;
const MANGLE_MASK = (1n << MANGLE_SHIFT) - 1n;
/** The seed reduction modulus `2^19937 - 20027` (`randmts.c:126`). */
const SEED_MODULUS = (1n << MANGLE_SHIFT) - 20027n;

/**
 * Reduce `r` modulo `2^19937 - 20023`, exactly as the `reduce:` label inside
 * `mangle_seed` does (`randmts.c:50-59`).
 *
 * Note this is *not* a canonical reduction: the loop stops as soon as
 * `r >> 19937 == 0`, so the result only satisfies `r < 2^19937` and may still
 * be `>= 2^19937 - 20023`.  Reproducing that quirk matters, because the top
 * bit of the result is what ends up in `mt[0]`.
 */
function mangleReduce(r: bigint): bigint {
  for (;;) {
    const t = r >> MANGLE_SHIFT;
    if (t === 0n) return r;
    r = (r & MANGLE_MASK) + t * 20023n;
  }
}

/**
 * Compute `r^1074888996 mod (2^19937 - 20023)`.
 *
 * Verbatim port of `mangle_seed` (GMP `rand/randmts.c:36-70`), including its
 * left-to-right square-and-multiply schedule (the leading bit `0x40000000` of
 * the exponent is consumed implicitly by starting from `r = b`).
 */
function mangleSeed(seed: bigint): bigint {
  const b = seed;
  let r = seed;
  let e = 0x40118124;
  let bit = 0x20000000;
  do {
    r = mangleReduce(r * r);
    if ((e & bit) !== 0) {
      e ^= bit;
      r = mangleReduce(r * b);
    }
    bit >>>= 1;
  } while (bit !== 0);
  return r;
}

/**
 * Random state container: a port of Sage's `randstate`.
 *
 * This is not a cryptographic RNG.  It matches Sage's centralized randstate
 * semantics — including the concrete bit stream for a given seed.
 */
export class RandState {
  /** The MT19937 state vector. */
  private readonly mt: Uint32Array;
  /** Index into the state vector; `MT_N` means "buffer exhausted". */
  private mti: number;
  /**
   * The seed this state was created (or reseeded) with.
   *
   * Sage keeps this separately from the generator state so that `seed()`
   * does not change as random numbers are drawn (`randstate.pyx:562-580`).
   */
  private _seed: bigint;
  /**
   * Lazily created CPython generator (`randstate.pyx:583-627`).
   */
  private _python_random: PythonRandom | null;

  constructor(seed?: bigint | number) {
    this.mt = new Uint32Array(MT_N);
    this.mti = MT_N;
    this._seed = 0n;
    this._python_random = null;
    this.set_seed(seed === undefined ? defaultSeed() : BigInt(seed));
  }

  /**
   * Set the RNG seed.
   *
   * Reference: `randstate.pyx:543-561`.  A seed of `0` is special: Sage skips
   * `gmp_randseed` entirely, so the generator keeps GMP's built-in default
   * buffer (`randmt.c:__gmp_randinit_mt_noseed`).
   */
  set_seed(seed: bigint | number): void {
    const seedValue = BigInt(seed);
    this._seed = seedValue;
    this._python_random = null;
    if (seedValue === 0n) {
      // gmp_randinit_default() state, untouched (randstate.pyx:556-559).
      for (let i = 0; i < MT_N; i++) {
        this.mt[i] = MT_DEFAULT_STATE[i]!;
      }
      this.mti = MT_WARM_UP % MT_N;
      return;
    }
    this._gmp_randseed(seedValue);
  }

  /**
   * GMP's `randseed_mt` (`rand/randmts.c:113-155`).
   */
  private _gmp_randseed(seed: bigint): void {
    // mpz_mod(seed1, seed, 2^19937 - 20027): always nonnegative.
    let seed1 = seed % SEED_MODULUS;
    if (seed1 < 0n) seed1 += SEED_MODULUS;
    seed1 += 2n; // "seed1 is now ready"
    seed1 = mangleSeed(seed1);

    // Copy the last bit into bit 31 of mt[0] and clear it.
    this.mt[0] = (seed1 >> 19936n) & 1n ? MT_UPPER_MASK : 0;
    seed1 &= (1n << 19936n) - 1n;

    // Split seed1 into N-1 32-bit chunks (mpz_export, least significant
    // first); `cnt` is the number of words written.
    let cnt = 1;
    let rest = seed1;
    while (rest > 0n) {
      this.mt[cnt++] = Number(rest & 0xffffffffn);
      rest >>= 32n;
    }
    while (cnt < MT_N) {
      this.mt[cnt++] = 0;
    }

    // Warm the generator up (WARM_UP = 2000, randmt.h:33).
    for (let i = 0; i < Math.floor(MT_WARM_UP / MT_N); i++) {
      mtRecalcBuffer(this.mt);
    }
    this.mti = MT_WARM_UP % MT_N;
  }

  /**
   * Return the initial seed of this random state.
   *
   * This is *not* the current state: it does not change when random numbers
   * are drawn.
   *
   * Reference: `randstate.pyx:562` `randstate.seed`.
   */
  seed(): bigint {
    return this._seed;
  }

  /**
   * GMP's `NEXT_RANDOM`: generate the next tempered 32-bit output word
   * (`rand/randmt.c:206-227`).
   */
  private next_uint32(): number {
    if (this.mti >= MT_N) {
      mtRecalcBuffer(this.mt);
      this.mti = 0;
    }
    return mtTemper(this.mt[this.mti++]!);
  }

  /**
   * Return a uniformly random integer in `[0, 2^nbits)`.
   *
   * This is GMP's `gmp_urandomb_ui`/`__gmp_randget_mt` (`rand/randmt.c:229-320`).
   * On a 64-bit build that consumes exactly `ceil(nbits/32)` generator
   * outputs, assembles them least-significant-word first and masks off the
   * excess high bits — the `GMP_NUMB_BITS == 64` branch splits the request
   * into `2*floor(nbits/64) + ceil((nbits mod 64)/32)` words, which is the
   * same count.
   */
  private urandomb(nbits: number): bigint {
    if (nbits === 0) {
      return 0n;
    }
    let result = 0n;
    let shift = 0n;
    for (let remaining = nbits; remaining > 0; remaining -= 32) {
      result |= BigInt(this.next_uint32()) << shift;
      shift += 32n;
    }
    return result & ((1n << BigInt(nbits)) - 1n);
  }

  /**
   * Return a 31-bit random number (0 <= x <= SAGE_RAND_MAX).
   *
   * Reference: `randstate.pyx:881` `randstate.c_random`.
   */
  c_random(): number {
    return Number(this.urandomb(31));
  }

  /**
   * Return a random floating-point number in [0, 1).
   *
   * Ported verbatim from `randstate.pyx:883-895`: it consumes *two* draws,
   * 25 bits scaled by 2^-25 plus 28 bits scaled by 2^-53.
   */
  c_rand_double(): number {
    const a = Number(this.urandomb(25)) * (1.0 / 33554432.0); // divide by 2^25
    const b = Number(this.urandomb(28)) * (1.0 / 9007199254740992.0); // divide by 2^53
    return a + b;
  }

  /**
   * Return a floating-point number in [0, 1).
   *
   * @deprecated Sage spells this `c_rand_double`; `random()` in
   * `sage.misc.randstate` is the module-level 31-bit integer generator. Kept
   * as an alias for existing callers.
   * @see Deviation: extra method not present in Sage.
   */
  random(): number {
    return this.c_rand_double();
  }

  /**
   * Return a uniformly random bigint in `[0, 2^bits)`.
   *
   * Mirrors GMP's `mpz_urandomb` (`mpz/urandomb.c`), which is just
   * `_gmp_rand` on `nbits` bits.
   */
  random_bits(bits: number): bigint {
    if (bits < 0) {
      throw new ValueError('number of bits must be nonnegative');
    }
    return this.urandomb(bits);
  }

  /**
   * Return a random `bits`-bit integer with **long runs** of equal bits.
   *
   * Port of GMP's `mpz_rrandomb` / `gmp_rrandomb` (`mpz/rrandomb.c`), which
   * `ZZ.random_element(x, distribution='mpz_rrandomb')` calls
   * (`integer_ring.pyx:822`).  This is a completely different distribution
   * from `mpz_urandomb`: the value starts as `2^bits - 1` and is then punched
   * with alternating runs, so a typical result sits within a few percent of
   * `2^bits`.
   *
   * ```c
   * i = BITS_TO_LIMBS (nbits) - 1;
   * rp[i] = GMP_NUMB_MAX >> (GMP_NUMB_BITS - (nbits % GMP_NUMB_BITS)) % GMP_NUMB_BITS;
   * for (i = i - 1; i >= 0; i--) rp[i] = GMP_NUMB_MAX;
   * _gmp_rand (&ranm, rstate, 32);
   * cap_chunksize = nbits / (ranm % 4 + 1);
   * cap_chunksize += cap_chunksize == 0;
   * bi = nbits;
   * for (;;) {
   *   _gmp_rand (&ranm, rstate, 32);
   *   chunksize = 1 + ranm % cap_chunksize;
   *   bi = (bi < chunksize) ? 0 : bi - chunksize;
   *   if (bi == 0) break;
   *   rp[bi / GMP_NUMB_BITS] ^= CNST_LIMB (1) << bi % GMP_NUMB_BITS;
   *   _gmp_rand (&ranm, rstate, 32);
   *   chunksize = 1 + ranm % cap_chunksize;
   *   bi = (bi < chunksize) ? 0 : bi - chunksize;
   *   mpn_incr_u (rp + bi / GMP_NUMB_BITS, CNST_LIMB (1) << bi % GMP_NUMB_BITS);
   *   if (bi == 0) break;
   * }
   * ```
   *
   * The limb layout is irrelevant here: XOR-ing bit `bi` and `mpn_incr_u` of
   * `1 << bi` are exactly `value ^= 2^bi` and `value += 2^bi` on the whole
   * integer, and the carry can never run past the bit cleared just before.
   */
  random_bits_rrandomb(bits: number): bigint {
    if (bits < 0) {
      throw new ValueError('number of bits must be nonnegative');
    }
    if (bits === 0) {
      return 0n;
    }

    let value = (1n << BigInt(bits)) - 1n;

    let ranm = this.urandomb(32);
    let capChunksize = Math.floor(bits / (Number(ranm % 4n) + 1));
    if (capChunksize === 0) {
      capChunksize = 1;
    }
    const cap = BigInt(capChunksize);

    let bi = bits;
    for (;;) {
      ranm = this.urandomb(32);
      let chunksize = 1 + Number(ranm % cap);
      bi = bi < chunksize ? 0 : bi - chunksize;
      if (bi === 0) {
        break;
      }

      value ^= 1n << BigInt(bi);

      ranm = this.urandomb(32);
      chunksize = 1 + Number(ranm % cap);
      bi = bi < chunksize ? 0 : bi - chunksize;

      value += 1n << BigInt(bi);

      if (bi === 0) {
        break;
      }
    }

    return value;
  }

  /**
   * A 128-bit integer suitable for seeding another random number generator.
   *
   * Port of `randstate.pyx:629-641`: `return ZZ.random_element(1<<128)`.
   *
   * Note the sharp edge Sage has here (shared with {@link python_random}):
   * `ZZ.random_element` reads `current_randstate()`, **not** `self`, so the
   * draw comes from the globally current state.  We reproduce that, including
   * the unconditional `den` draw `ZZ.random_element` burns first
   * (`integer_ring.pyx:801`).
   */
  ZZ_seed(): bigint {
    const source = current_randstate();
    source.c_random(); // `integer_ring.pyx:801` `den`, computed and discarded
    return source.random_below(1n << 128n);
  }

  /**
   * The same 128-bit seed as {@link ZZ_seed}.
   *
   * Port of `randstate.pyx:643-656`: `return int(ZZ.random_element(1<<128))`.
   * Python's `int` and Sage's `Integer` are both plain bigints for us, so the
   * two methods coincide; both are kept because Sage exposes both.
   */
  long_seed(): bigint {
    return this.ZZ_seed();
  }

  /**
   * Return a random bigint in [0, n).
   *
   * Port of `mpz_urandomm` (GMP `mpz/urandomm.c:33-95`): draw
   * `bitlen(n) - [n is a power of two]` bits and reject while the value is
   * `>= n`, giving up after `MAX_URANDOMM_ITER = 80` tries and subtracting
   * `n` once.  Sage reaches this through `ZZ.random_element(n)`
   * (`integer_ring.pyx:811`).
   */
  random_below(n: bigint): bigint {
    if (n <= 0n) {
      throw new ValueError('n must be positive');
    }
    if (n === 1n) {
      return 0n;
    }
    const isPow2 = (n & (n - 1n)) === 0n;
    const nbits = n.toString(2).length - (isPow2 ? 1 : 0);
    let candidate = 0n;
    for (let count = MAX_URANDOMM_ITER; count > 0; count--) {
      candidate = this.urandomb(nbits);
      if (candidate < n) {
        return candidate;
      }
    }
    // Too many iterations; return result mod n == result - n.
    return candidate - n;
  }

  /**
   * Return a random bigint in [min, max] inclusive.
   */
  randint(min: bigint, max: bigint): bigint {
    if (min > max) {
      throw new ValueError('min must be <= max');
    }
    const range = max - min + 1n;
    return min + this.random_below(range);
  }

  /**
   * Return this state's CPython `random.Random` object.
   *
   * Port of `randstate.python_random` (`randstate.pyx:583-627`): the first
   * call creates a `random.Random` seeded with `int(ZZ.random_element(1<<128))`.
   * `ZZ.random_element(x)` first burns one `c_random()` draw on its unused
   * `den` variable (`integer_ring.pyx:801`) and then calls `mpz_urandomm`
   * (`integer_ring.pyx:811`).
   *
   * @see Deviation: none — but note the sharp edge Sage has here.
   * `ZZ.random_element` reads `current_randstate()`, *not* `self`, so the
   * derived seed comes from the globally current state rather than from this
   * one.  We reproduce that: `set_random_seed(0)` followed by
   * `new RandState(314159).python_random().random()` gives Sage's
   * `0.111439293741037` (`randstate.pyx:578`), which is a seed-0 value.  Pass
   * an explicit `seed` when you want a self-contained stream.
   *
   * @param seed - explicit seed, as Sage's `seed=` keyword
   */
  python_random(seed?: bigint): PythonRandom {
    // `randstate.pyx:617` returns the cached generator BEFORE `seed` is ever
    // looked at, so once one exists an explicit `seed=` is silently discarded
    // and the stream simply continues.  The cache check must therefore be
    // unconditional.
    if (this._python_random !== null) {
      return this._python_random;
    }
    const rand = new PythonRandom();
    if (seed === undefined) {
      // ZZ.random_element(1 << 128), evaluated on the *current* randstate.
      const source = current_randstate();
      source.c_random(); // integer_ring.pyx:801 `den`, discarded
      rand.seed(source.random_below(1n << 128n));
    } else {
      rand.seed(seed);
    }
    this._python_random = rand;
    return rand;
  }
}

/**
 * CPython's `random.Random`.
 *
 * Port of `Modules/_randommodule.c` (`init_genrand`, `init_by_array`,
 * `genrand_uint32`, `random_seed`, `random_random`, `getrandbits`) and of the
 * pure-Python layer in `Lib/random.py` (`_randbelow_with_getrandbits`,
 * `randrange`, `randint`, `normalvariate`).
 *
 * Sage exposes this generator as `randstate.python_random()`
 * (`randstate.pyx:583`) and routes `sage.misc.prandom` through it, so it is
 * the generator behind `IntegerModRing.random_element()`
 * (`integer_mod_ring.py:1543-1547`) and `DiscreteGaussianDistributionLatticeSampler`'s
 * `normalvariate` calls (`discrete_gaussian_lattice.py:891`).
 */
export class PythonRandom {
  private readonly mt: Uint32Array;
  private mti: number;

  constructor(seed?: bigint | number) {
    this.mt = new Uint32Array(MT_N);
    this.mti = MT_N + 1;
    this.seed(seed === undefined ? defaultSeed() : BigInt(seed));
  }

  /**
   * `random_seed` (`_randommodule.c`): use `abs(n)` split into 32-bit words,
   * least significant first, and run `init_by_array`.
   */
  seed(n: bigint | number): void {
    let s = BigInt(n);
    if (s < 0n) s = -s;
    const key: number[] = [];
    while (s > 0n) {
      key.push(Number(s & 0xffffffffn));
      s >>= 32n;
    }
    if (key.length === 0) key.push(0);
    this._init_by_array(key);
  }

  /** MT19937 `init_genrand`. */
  private _init_genrand(s: number): void {
    this.mt[0] = s >>> 0;
    for (let i = 1; i < MT_N; i++) {
      const prev = this.mt[i - 1]!;
      this.mt[i] = (Math.imul(1812433253, prev ^ (prev >>> 30)) + i) >>> 0;
    }
    this.mti = MT_N;
  }

  /** MT19937 `init_by_array`. */
  private _init_by_array(key: number[]): void {
    this._init_genrand(19650218);
    let i = 1;
    let j = 0;
    let k = Math.max(MT_N, key.length);
    for (; k > 0; k--) {
      const prev = this.mt[i - 1]!;
      this.mt[i] =
        (((this.mt[i]! ^ Math.imul(prev ^ (prev >>> 30), 1664525)) >>> 0) + key[j]! + j) >>> 0;
      i++;
      j++;
      if (i >= MT_N) {
        this.mt[0] = this.mt[MT_N - 1]!;
        i = 1;
      }
      if (j >= key.length) {
        j = 0;
      }
    }
    for (k = MT_N - 1; k > 0; k--) {
      const prev = this.mt[i - 1]!;
      this.mt[i] = (((this.mt[i]! ^ Math.imul(prev ^ (prev >>> 30), 1566083941)) >>> 0) - i) >>> 0;
      i++;
      if (i >= MT_N) {
        this.mt[0] = this.mt[MT_N - 1]!;
        i = 1;
      }
    }
    // MSB is 1, assuring a non-zero initial array.
    this.mt[0] = MT_UPPER_MASK;
    this.mti = MT_N;
  }

  /** `genrand_uint32` (`_randommodule.c`). */
  private genrand_uint32(): number {
    if (this.mti >= MT_N) {
      mtRecalcBuffer(this.mt);
      this.mti = 0;
    }
    return mtTemper(this.mt[this.mti++]!);
  }

  /**
   * `random_random`: `(a*67108864.0 + b) * 2^-53` with `a` the top 27 bits and
   * `b` the top 26 bits of two consecutive outputs.
   */
  random(): number {
    const a = this.genrand_uint32() >>> 5;
    const b = this.genrand_uint32() >>> 6;
    return (a * 67108864.0 + b) * (1.0 / 9007199254740992.0);
  }

  /**
   * `getrandbits(k)`: `ceil(k/32)` outputs, least significant word first, with
   * the *most* significant word right-shifted to drop the surplus low bits.
   */
  getrandbits(k: number): bigint {
    if (k < 0) {
      throw new ValueError('number of bits must be non-negative');
    }
    if (k === 0) {
      return 0n;
    }
    if (k <= 32) {
      return BigInt(this.genrand_uint32() >>> (32 - k));
    }
    const words = Math.floor((k - 1) / 32) + 1;
    let result = 0n;
    let remaining = k;
    for (let i = 0; i < words; i++, remaining -= 32) {
      let r = this.genrand_uint32();
      if (remaining < 32) {
        r >>>= 32 - remaining;
      }
      result |= BigInt(r >>> 0) << BigInt(32 * i);
    }
    return result;
  }

  /**
   * `Random._randbelow_with_getrandbits` (`Lib/random.py`).
   */
  private _randbelow(n: bigint): bigint {
    if (n <= 0n) {
      return 0n;
    }
    const k = n.toString(2).length;
    let r = this.getrandbits(k);
    while (r >= n) {
      r = this.getrandbits(k);
    }
    return r;
  }

  /**
   * `Random.randrange(start, stop)` for the two-argument, unit-step form.
   */
  randrange(start: bigint | number, stop?: bigint | number): bigint {
    const a = BigInt(start);
    if (stop === undefined) {
      if (a <= 0n) {
        // CPython 3.11 `Lib/random.py`: `"empty range for randrange(%d)"` in
        // 3.10 and earlier, bare `"empty range for randrange()"` from 3.11.
        // SageMath 10.3 runs on 3.11.8.
        throw new ValueError('empty range for randrange()');
      }
      return this._randbelow(a);
    }
    const b = BigInt(stop);
    const width = b - a;
    if (width <= 0n) {
      // CPython 3.11 `Lib/random.py`:
      // `"empty range for randrange() (%d, %d, %d)" % (istart, istop, width)`.
      throw new ValueError(`empty range for randrange() (${a}, ${b}, ${width})`);
    }
    return a + this._randbelow(width);
  }

  /**
   * `Random.randint(a, b)`: `randrange(a, b+1)`.
   */
  randint(a: bigint | number, b: bigint | number): bigint {
    return this.randrange(BigInt(a), BigInt(b) + 1n);
  }

  /**
   * `Random.normalvariate(mu, sigma)` (`Lib/random.py`): Kinderman and Monahan
   * ratio-of-uniforms method, with `NV_MAGICCONST = 4*exp(-0.5)/sqrt(2)`.
   */
  normalvariate(mu = 0.0, sigma = 1.0): number {
    const NV_MAGICCONST = (4 * Math.exp(-0.5)) / Math.sqrt(2.0);
    for (;;) {
      const u1 = this.random();
      const u2 = 1.0 - this.random();
      const z = (NV_MAGICCONST * (u1 - 0.5)) / u2;
      const zz = (z * z) / 4.0;
      if (zz <= -Math.log(u2)) {
        return mu + z * sigma;
      }
    }
  }
}

/**
 * Return the current global random state.
 */
export function current_randstate(): RandState {
  if (_currentRandState === null) {
    _currentRandState = new RandState();
  }
  return _currentRandState;
}

/**
 * Set the global random seed.
 *
 * Reference: `randstate.pyx:978` `set_random_seed`.
 */
export function set_random_seed(seed?: number | bigint): void {
  if (_currentRandState === null) {
    _currentRandState = new RandState(seed);
  } else if (seed === undefined) {
    _currentRandState.set_seed(defaultSeed());
  } else {
    _currentRandState.set_seed(seed);
  }
}

/**
 * Return a 31-bit random number.
 *
 * Intended as a drop-in replacement for the libc `random()` function.
 *
 * Reference: `randstate.pyx:1008` `random`.
 */
export function random(): bigint {
  return BigInt(current_randstate().c_random());
}

/**
 * Return the initial seed used to create the current random state.
 *
 * Reference: `randstate.pyx:1024` `initial_seed`.
 */
export function initial_seed(): bigint {
  return current_randstate().seed();
}
