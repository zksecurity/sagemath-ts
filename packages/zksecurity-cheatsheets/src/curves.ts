import { GF } from '../../sagemath-ts/src/rings/finite_rings/index.js';
import { EllipticCurve } from '../../sagemath-ts/src/schemes/elliptic_curves/index.js';

// =============================================================================
// Interfaces
// =============================================================================

export interface PrimeFieldCurveParams {
  name: string;
  field: {
    p: bigint;
  };
  a: bigint;
  b: bigint;
  order: bigint;
  cofactor?: bigint;
  generator: {
    x: bigint;
    y: bigint;
  };
}

export interface MontgomeryCurveParams {
  name: string;
  field: {
    p: bigint;
  };
  A: bigint;
  B: bigint;
  order: bigint;
  cofactor: bigint;
  generator: {
    u: bigint;
    v: bigint;
  };
}

export interface TwistedEdwardsCurveParams {
  name: string;
  field: {
    p: bigint;
  };
  a: bigint;
  d: bigint;
  order: bigint;
  cofactor: bigint;
  subgroupOrder: bigint;
  generator: {
    x: bigint;
    y: bigint;
  };
}

export interface Bls12_381Params {
  name: string;
  field: {
    p: bigint;
  };
  r: bigint;
  a: bigint;
  b: bigint;
  g1: {
    x: bigint;
    y: bigint;
    cofactor?: bigint;
  };
  g2: {
    x: [bigint, bigint];
    y: [bigint, bigint];
    cofactor?: bigint;
  };
}

export interface Bn254Params {
  name: string;
  field: {
    p: bigint;
  };
  r: bigint;
  a: bigint;
  b: bigint;
  g1: {
    x: bigint;
    y: bigint;
    cofactor?: bigint;
  };
  g2: {
    x: [bigint, bigint];
    y: [bigint, bigint];
    cofactor?: bigint;
  };
}

// =============================================================================
// Utilities
// =============================================================================

function hexToBigInt(hex: string): bigint {
  return BigInt(`0x${hex.replace(/\s+/g, '')}`);
}

export function createPrimeFieldCurve(params: PrimeFieldCurveParams) {
  const field = GF(params.field.p, { check: false });
  const curve = EllipticCurve(field, [params.a, params.b]);
  const generator = curve.point(params.generator.x, params.generator.y);

  return { field, curve, generator };
}

// =============================================================================
// Standard Curves (Weierstrass form)
// =============================================================================

export const secp256k1: PrimeFieldCurveParams = {
  name: 'secp256k1',
  field: {
    p: hexToBigInt('FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F'),
  },
  a: 0n,
  b: 7n,
  order: hexToBigInt('FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141'),
  cofactor: 1n,
  generator: {
    x: hexToBigInt('79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798'),
    y: hexToBigInt('483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8'),
  },
};

export const p256: PrimeFieldCurveParams = {
  name: 'secp256r1',
  field: {
    p: hexToBigInt('FFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFF'),
  },
  a: hexToBigInt('FFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFC'),
  b: hexToBigInt('5AC635D8AA3A93E7B3EBBD55769886BC651D06B0CC53B0F63BCE3C3E27D2604B'),
  order: hexToBigInt('FFFFFFFF00000000FFFFFFFFFFFFFFFFBCE6FAADA7179E84F3B9CAC2FC632551'),
  cofactor: 1n,
  generator: {
    x: hexToBigInt('6B17D1F2E12C4247F8BCE6E563A440F277037D812DEB33A0F4A13945D898C296'),
    y: hexToBigInt('4FE342E2FE1A7F9B8EE7EB4A7C0F9E162BCE33576B315ECECBB6406837BF51F5'),
  },
};

// =============================================================================
// Pairing Curves
// =============================================================================

export const bls12_381: Bls12_381Params = {
  name: 'bls12-381',
  field: {
    p: hexToBigInt(
      '1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f624' +
        '1eabfffeb153ffffb9feffffffffaaab'
    ),
  },
  r: hexToBigInt('73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001'),
  a: 0n,
  b: 4n,
  g1: {
    x: BigInt(
      '3685416753713387016781088315183077757961620795782546409894578378688607592378376318836054947676345821548104185464507'
    ),
    y: BigInt(
      '1339506544944476473020471379941921221584933875938349620426543736416511423956333506472724655353366534992391756441569'
    ),
  },
  g2: {
    x: [
      BigInt(
        '352701069587466618187139116011060144890029952792775240219908644239793785735715026873347600343865175952761926303160'
      ),
      BigInt(
        '3059144344244213709971259814753781636986470325476647558659373206291635324768958432433509563104347017837885763365758'
      ),
    ],
    y: [
      BigInt(
        '1985150602287291935568054521177171638300868978215655730859378665066344726373823718423869104263333984641494340347905'
      ),
      BigInt(
        '927553665492332455747201965776037880757740193453592970025027978793976877002675564980949289727957565575433344219582'
      ),
    ],
  },
};

export const bls12_381_g1: PrimeFieldCurveParams = {
  name: 'bls12-381-g1',
  field: {
    p: bls12_381.field.p,
  },
  a: bls12_381.a,
  b: bls12_381.b,
  order: bls12_381.r,
  generator: {
    x: bls12_381.g1.x,
    y: bls12_381.g1.y,
  },
};

/**
 * BN254 (alt_bn128) - Ethereum's pairing curve for precompiles.
 * Used by Groth16 verifiers on-chain (ecAdd, ecMul, ecPairing).
 * Curve: y² = x³ + 3
 */
export const bn254: Bn254Params = {
  name: 'bn254',
  field: {
    p: 21888242871839275222246405745257275088696311157297823662689037894645226208583n,
  },
  r: 21888242871839275222246405745257275088548364400416034343698204186575808495617n,
  a: 0n,
  b: 3n,
  g1: {
    x: 1n,
    y: 2n,
  },
  g2: {
    // Fq2 tower with u² + 1 = 0
    x: [
      10857046999023057135944570762232829481370756359578518086990519993285655852781n,
      11559732032986387107991004021392285783925812861821192530917403151452391805634n,
    ],
    y: [
      8495653923123431417604973247489272438418190587263600148770280649306958101930n,
      4082367875863433681332203403145435568316851327593401208105741076214120093531n,
    ],
  },
};

export const bn254_g1: PrimeFieldCurveParams = {
  name: 'bn254-g1',
  field: {
    p: bn254.field.p,
  },
  a: bn254.a,
  b: bn254.b,
  order: bn254.r,
  generator: {
    x: bn254.g1.x,
    y: bn254.g1.y,
  },
};

// =============================================================================
// Curve Cycles (for recursive SNARKs)
// =============================================================================

/**
 * Grumpkin - companion curve to BN254, forming a 2-cycle.
 * Used for recursive SNARKs (Aztec/Noir).
 * Curve: y² = x³ - 1
 *
 * Cycle relationship:
 * - Grumpkin.p = BN254.r (BN254's scalar field is Grumpkin's base field)
 * - Grumpkin.order = BN254.p (BN254's base field is Grumpkin's scalar field)
 */
export const grumpkin: PrimeFieldCurveParams = {
  name: 'grumpkin',
  field: {
    p: 21888242871839275222246405745257275088548364400416034343698204186575808495617n,
  },
  a: 0n,
  b: 21888242871839275222246405745257275088548364400416034343698204186575808495616n, // -1 mod p
  order: 21888242871839275222246405745257275088696311157297823662689037894645226208583n,
  cofactor: 1n,
  generator: {
    x: 1n,
    y: 17631683881184975370165255887551781615748388533673675138860n,
  },
};

/**
 * Pasta curves - a 2-cycle used by Halo2, Mina, and Zcash.
 * Pallas.order = Vesta.p and Vesta.order = Pallas.p
 */
export const pallas: PrimeFieldCurveParams = {
  name: 'pallas',
  field: {
    p: 0x40000000000000000000000000000000224698fc094cf91b992d30ed00000001n,
  },
  a: 0n,
  b: 5n,
  order: 0x40000000000000000000000000000000224698fc0994a8dd8c46eb2100000001n,
  cofactor: 1n,
  generator: {
    x: 0x40000000000000000000000000000000224698fc094cf91b992d30ed00000001n - 1n, // -1 mod p
    y: 2n,
  },
};

export const vesta: PrimeFieldCurveParams = {
  name: 'vesta',
  field: {
    p: 0x40000000000000000000000000000000224698fc0994a8dd8c46eb2100000001n,
  },
  a: 0n,
  b: 5n,
  order: 0x40000000000000000000000000000000224698fc094cf91b992d30ed00000001n,
  cofactor: 1n,
  generator: {
    x: 0x40000000000000000000000000000000224698fc0994a8dd8c46eb2100000001n - 1n, // -1 mod p
    y: 2n,
  },
};

// =============================================================================
// Montgomery Curves
// =============================================================================

const curve25519_p = 2n ** 255n - 19n;
const curve25519_l = 7237005577332262213973186563042994240857116359379907606001950938285454250989n;

/**
 * Curve25519 - Montgomery curve for X25519 key exchange.
 * Equation: v² = u³ + 486662u² + u
 */
export const curve25519: MontgomeryCurveParams = {
  name: 'curve25519',
  field: {
    p: curve25519_p,
  },
  A: 486662n,
  B: 1n,
  order: 8n * curve25519_l,
  cofactor: 8n,
  generator: {
    u: 9n,
    v: 14781619447589544791020593568409986887264606134616475288964881837755586237401n,
  },
};

// =============================================================================
// Twisted Edwards Curves
// =============================================================================

/**
 * Ed25519 - Twisted Edwards curve for EdDSA signatures.
 * Equation: -x² + y² = 1 + dx²y²
 * Birationally equivalent to Curve25519.
 */
export const ed25519: TwistedEdwardsCurveParams = {
  name: 'ed25519',
  field: {
    p: curve25519_p,
  },
  a: curve25519_p - 1n, // -1 mod p
  d: 37095705934669439343138083508754565189542113879843219016388785533085940283555n,
  order: 8n * curve25519_l,
  cofactor: 8n,
  subgroupOrder: curve25519_l,
  generator: {
    x: 15112221349535807912866137220509078935008241517919376607934839071348814711696n,
    y: 46316835694926478169428394003475163141307993866256225615783033603165251855960n,
  },
};

/**
 * BabyJubJub - Twisted Edwards curve embedded over BN254's scalar field.
 * Used for in-circuit EdDSA signatures and Pedersen hashes in Ethereum ZK systems.
 * Equation: 168700x² + y² = 1 + 168696x²y²
 */
export const babyJubJub: TwistedEdwardsCurveParams = {
  name: 'babyjubjub',
  field: {
    p: 21888242871839275222246405745257275088548364400416034343698204186575808495617n, // BN254.r
  },
  a: 168700n,
  d: 168696n,
  order: 21888242871839275222246405745257275088614511777268538073601725287587578984328n,
  cofactor: 8n,
  subgroupOrder: 2736030358979909402780800718157159386076813972158567259200215660948447373041n,
  generator: {
    x: 5299619240641551281634865583518297030282874472190772894086521144482721001553n,
    y: 16950150798460657717958625567821834550301663161624707787222815936182638968203n,
  },
};

// =============================================================================
// Exports
// =============================================================================

export const curves = {
  // Standard
  secp256k1,
  p256,
  // Pairing
  bls12_381,
  bls12_381_g1,
  bn254,
  bn254_g1,
  // Cycles
  grumpkin,
  pallas,
  vesta,
  // Montgomery
  curve25519,
  // Edwards
  ed25519,
  babyJubJub,
};
