/** Port of the original independent fixture builder. Uses native RNGs and affine EC,
 * not the timestamp engine or its Jacobian curve implementation. Generated output is ephemeral.
 */
import { createHash } from 'node:crypto';
import { json } from '../../../src/timestampseeds/json.js';
import { References } from './references.js';

const CURVES = {
  p256: [
    0xffffffff00000001000000000000000000000000ffffffffffffffffffffffffn,
    -3n,
    0xffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc632551n,
    [
      0x6b17d1f2e12c4247f8bce6e563a440f277037d812deb33a0f4a13945d898c296n,
      0x4fe342e2fe1a7f9b8ee7eb4a7c0f9e162bce33576b315ececbb6406837bf51f5n,
    ],
    '3059301306072a8648ce3d020106082a8648ce3d03010703420004',
  ],
  secp256k1: [
    (1n << 256n) - (1n << 32n) - 977n,
    0n,
    0xfffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141n,
    [
      0x79be667ef9dcbbac55a06295ce870b07029bfcdb2dce28d959f2815b16f81798n,
      0x483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8n,
    ],
    '3056301006072a8648ce3d020106052b8104000a03420004',
  ],
} as const;
const RNG_NAMES = [
  'lrand48',
  'drand48-state',
  'glibc-rand',
  'ansi-lcg31',
  'python-mt19937',
  'java-random',
  'msvc-rand',
  'cpp-mt19937',
];
type Point = readonly [bigint, bigint] | null;
function multiply(value: bigint, prime: bigint, a: bigint, point: Point): Point {
  const mod = (n: bigint) => ((n % prime) + prime) % prime;
  function inverse(n: bigint): bigint {
    let b = prime,
      x = 1n,
      y = 0n;
    n = mod(n);
    while (b) {
      const q = n / b;
      [n, b] = [b, n - q * b];
      [x, y] = [y, x - q * y];
    }
    if (n !== 1n) throw new Error('Noninvertible reference coordinate');
    return mod(x);
  }
  function add(left: Point, right: Point): Point {
    if (!left) return right;
    if (!right) return left;
    const [x, y] = left,
      [u, v] = right;
    if (x === u && mod(y + v) === 0n) return null;
    const slope = mod(
      x === u && y === v ? (3n * x * x + a) * inverse(2n * y) : (v - y) * inverse(u - x)
    );
    const rx = mod(slope * slope - x - u);
    return [rx, mod(slope * (x - rx) - y)];
  }
  let result: Point = null;
  while (value) {
    if (value & 1n) result = add(result, point);
    point = add(point, point);
    value >>= 1n;
  }
  return result;
}
export function public_reference(curve: keyof typeof CURVES, value: bigint, profile: string) {
  const [prime, a, order, generator, prefix] = CURVES[curve];
  const scalar = profile === 'concat-be-mod' ? (value % (order - 1n)) + 1n : value;
  if (scalar < 1n || scalar >= order) return null;
  const [x, y] = multiply(scalar, prime, a, generator)!;
  const hexX = x.toString(16).padStart(64, '0'),
    hexY = y.toString(16).padStart(64, '0');
  return {
    sec1_compressed_hex: (y & 1n ? '03' : '02') + hexX,
    spki_sha256: createHash('sha256')
      .update(Buffer.from(prefix + hexX + hexY, 'hex'))
      .digest('hex'),
  };
}
export function build(reference: References) {
  const cases: [string, bigint, string, number, number][] = RNG_NAMES.map((rng) => [
    rng,
    42n,
    'concat-be-reject',
    0,
    0,
  ]);
  cases.push(
    ...RNG_NAMES.map((rng): [string, bigint, string, number, number] => [
      rng,
      1700000000n,
      'low-byte-reject',
      3,
      2,
    ]),
    ['drand48-state', 0xffffffffn, 'concat-le-reject', 7, 1],
    ['cpp-mt19937', 42n, 'buffer-le-reject', 1, 1],
    ['cpp-mt19937', 2n, 'concat-be-mod', 0, 0],
    ['python-mt19937', (1n << 64n) + 42n, 'python-getrandbits-reject', 5, 3],
    ['python-mt19937', (1n << 64n) + 42n, 'python-randbytes-reject', 5, 3],
    ['java-random', (1n << 48n) + 42n, 'java-nextbytes-reject', 1, 2]
  );
  return {
    reference: 'tests/reference/build_fixtures.ts; native references plus independent affine EC',
    cases: cases.flatMap(([rng, seed, profile, skip_draws, key_index]) => {
      const value = reference.value(rng, seed, profile, skip_draws, key_index);
      return (['p256', 'secp256k1'] as const).map((curve) => ({
        rng,
        seed,
        profile,
        curve,
        skip_draws,
        key_index,
        public: public_reference(curve, value, profile),
      }));
    }),
  };
}
if (import.meta.main) {
  const reference = new References();
  try {
    console.log(json(build(reference)));
  } finally {
    reference.close();
  }
}
