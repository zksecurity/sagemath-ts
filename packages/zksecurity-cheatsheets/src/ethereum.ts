export interface PrecompileInfo {
  address: string;
  name: string;
  description: string;
  gasCost: string | number; // Can be formula or fixed
  inputFormat: string;
  outputFormat: string;
}

export const precompiles: Record<string, PrecompileInfo> = {
  ecRecover: {
    address: '0x01',
    name: 'ecRecover',
    description: 'Recover signer address from secp256k1 signature',
    gasCost: 3000,
    inputFormat: 'hash (32) + v (32) + r (32) + s (32)',
    outputFormat: 'address (32, left-padded)',
  },
  sha256: {
    address: '0x02',
    name: 'SHA256',
    description: 'SHA-256 hash',
    gasCost: '60 + 12 * ceil(len/32)',
    inputFormat: 'arbitrary bytes',
    outputFormat: 'hash (32)',
  },
  ripemd160: {
    address: '0x03',
    name: 'RIPEMD-160',
    description: 'RIPEMD-160 hash',
    gasCost: '600 + 120 * ceil(len/32)',
    inputFormat: 'arbitrary bytes',
    outputFormat: 'hash (32, left-padded)',
  },
  identity: {
    address: '0x04',
    name: 'Identity',
    description: 'Data copy',
    gasCost: '15 + 3 * ceil(len/32)',
    inputFormat: 'arbitrary bytes',
    outputFormat: 'same bytes',
  },
  modexp: {
    address: '0x05',
    name: 'ModExp',
    description: 'Modular exponentiation',
    gasCost: 'complex formula based on inputs',
    inputFormat: 'Bsize (32) + Esize (32) + Msize (32) + B + E + M',
    outputFormat: 'result (Msize bytes)',
  },
  ecAdd: {
    address: '0x06',
    name: 'BN254 ecAdd',
    description: 'BN254 G1 point addition',
    gasCost: 150, // After Berlin
    inputFormat: 'P1.x (32) + P1.y (32) + P2.x (32) + P2.y (32)',
    outputFormat: 'P3.x (32) + P3.y (32)',
  },
  ecMul: {
    address: '0x07',
    name: 'BN254 ecMul',
    description: 'BN254 G1 scalar multiplication',
    gasCost: 6000, // After Berlin
    inputFormat: 'P.x (32) + P.y (32) + scalar (32)',
    outputFormat: 'Q.x (32) + Q.y (32)',
  },
  ecPairing: {
    address: '0x08',
    name: 'BN254 ecPairing',
    description: 'BN254 pairing check',
    gasCost: '45000 + 34000 * k', // k = number of pairs
    inputFormat:
      'k * (G1.x (32) + G1.y (32) + G2.x[0] (32) + G2.x[1] (32) + G2.y[0] (32) + G2.y[1] (32))',
    outputFormat: '0x01 if pairing check passes, else 0x00 (32 bytes)',
  },
  blake2f: {
    address: '0x09',
    name: 'BLAKE2b F',
    description: 'BLAKE2b compression function',
    gasCost: '1 * rounds',
    inputFormat: 'rounds (4) + h (64) + m (128) + t (16) + f (1)',
    outputFormat: 'h (64)',
  },
  pointEvaluation: {
    address: '0x0a',
    name: 'Point Evaluation (EIP-4844)',
    description: 'KZG point evaluation for blob transactions',
    gasCost: 50000,
    inputFormat:
      'versioned_hash (32) + z (32) + y (32) + commitment (48) + proof (48)',
    outputFormat: 'FIELD_ELEMENTS_PER_BLOB (32) + BLS_MODULUS (32)',
  },
};

// BN254 point encoding helpers
export const bn254Encoding = {
  // Infinity point
  G1_INFINITY: '0x' + '00'.repeat(64),
  // Generator (1, 2)
  G1_GENERATOR:
    '0x' +
    '0000000000000000000000000000000000000000000000000000000000000001' +
    '0000000000000000000000000000000000000000000000000000000000000002',
  // Field modulus
  FIELD_MODULUS:
    21888242871839275222246405745257275088696311157297823662689037894645226208583n,
  // Scalar field modulus
  SCALAR_MODULUS:
    21888242871839275222246405745257275088548364400416034343698204186575808495617n,
};

export const ethereum = {
  precompiles,
  bn254Encoding,
};
