export interface PrimeFieldParams {
  name: string;
  p: bigint;
  bits: number;
  twoAdicity: number; // largest k where 2^k divides p-1
  generator: bigint; // multiplicative generator
  description: string;
}

// Goldilocks field - used by Plonky2, Polygon Zero
// p = 2^64 - 2^32 + 1
export const goldilocks: PrimeFieldParams = {
  name: 'goldilocks',
  p: 2n ** 64n - 2n ** 32n + 1n, // 18446744069414584321
  bits: 64,
  twoAdicity: 32,
  generator: 7n,
  description: 'Plonky2, Polygon Zero - efficient 64-bit arithmetic',
};

// Mersenne31 - used by Plonky3, Circle STARKs
// p = 2^31 - 1
export const mersenne31: PrimeFieldParams = {
  name: 'mersenne31',
  p: 2n ** 31n - 1n, // 2147483647
  bits: 31,
  twoAdicity: 1,
  generator: 7n,
  description: 'Plonky3, Circle STARKs - fast modular reduction',
};

// BabyBear - used by RISC Zero, Plonky3
// p = 2^31 - 2^27 + 1
export const babyBear: PrimeFieldParams = {
  name: 'babybear',
  p: 2n ** 31n - 2n ** 27n + 1n, // 2013265921
  bits: 31,
  twoAdicity: 27,
  generator: 31n,
  description: 'RISC Zero, Plonky3 - high 2-adicity for FFTs',
};

// Koala Bear - used by some newer systems
// p = 2^31 - 2^24 + 1
export const koalaBear: PrimeFieldParams = {
  name: 'koalabear',
  p: 2n ** 31n - 2n ** 24n + 1n, // 2130706433
  bits: 31,
  twoAdicity: 24,
  generator: 3n,
  description: 'Alternative to BabyBear with different 2-adicity',
};

// BN254 scalar field
export const bn254Scalar: PrimeFieldParams = {
  name: 'bn254-scalar',
  p: 21888242871839275222246405745257275088548364400416034343698204186575808495617n,
  bits: 254,
  twoAdicity: 28,
  generator: 5n,
  description: 'BN254 scalar field (Fr) - Ethereum ZK circuits',
};

// BLS12-381 scalar field
export const bls12_381Scalar: PrimeFieldParams = {
  name: 'bls12-381-scalar',
  p: 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001n,
  bits: 255,
  twoAdicity: 32,
  generator: 7n,
  description: 'BLS12-381 scalar field (Fr) - Zcash, Ethereum 2.0',
};

// Pasta scalar fields (they are each other's base fields)
export const pallasScalar: PrimeFieldParams = {
  name: 'pallas-scalar',
  p: 0x40000000000000000000000000000000224698fc0994a8dd8c46eb2100000001n,
  bits: 255,
  twoAdicity: 32,
  generator: 5n,
  description: 'Pallas scalar field = Vesta base field',
};

export const vestaScalar: PrimeFieldParams = {
  name: 'vesta-scalar',
  p: 0x40000000000000000000000000000000224698fc094cf91b992d30ed00000001n,
  bits: 255,
  twoAdicity: 32,
  generator: 5n,
  description: 'Vesta scalar field = Pallas base field',
};

export const fields = {
  goldilocks,
  mersenne31,
  babyBear,
  koalaBear,
  bn254Scalar,
  bls12_381Scalar,
  pallasScalar,
  vestaScalar,
};
