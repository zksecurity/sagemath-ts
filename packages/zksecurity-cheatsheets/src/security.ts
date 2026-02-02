export interface SecurityLevel {
  curve: string;
  bits: number; // Classical security bits
  quantumBits?: number; // Post-quantum security (if applicable)
  keySize: number; // Private key size in bytes
  signatureSize?: number; // Signature size in bytes (if applicable)
  notes: string;
}

export const securityLevels: Record<string, SecurityLevel> = {
  // ECDSA/Schnorr curves
  secp256k1: {
    curve: 'secp256k1',
    bits: 128,
    keySize: 32,
    signatureSize: 64,
    notes: 'Bitcoin, Ethereum L1',
  },
  p256: {
    curve: 'secp256r1',
    bits: 128,
    keySize: 32,
    signatureSize: 64,
    notes: 'NIST standard, WebAuthn, TLS',
  },
  ed25519: {
    curve: 'ed25519',
    bits: 128,
    keySize: 32,
    signatureSize: 64,
    notes: 'Solana, Sui, fast signatures',
  },

  // Pairing curves
  bn254: {
    curve: 'bn254',
    bits: 100, // Reduced due to NFS attacks on small embedding degree
    keySize: 32,
    notes: 'Ethereum precompiles - security margin decreased',
  },
  bls12_381: {
    curve: 'bls12-381',
    bits: 128,
    keySize: 32,
    signatureSize: 48, // Compressed G1
    notes: 'Zcash, Ethereum 2.0, recommended pairing curve',
  },

  // ZK circuit curves
  pallas: {
    curve: 'pallas',
    bits: 128,
    keySize: 32,
    notes: 'Halo2, Mina - part of Pasta cycle',
  },
  grumpkin: {
    curve: 'grumpkin',
    bits: 128,
    keySize: 32,
    notes: 'Aztec/Noir - companion to BN254',
  },
};

// Recommended curves by use case
export const recommendations = {
  signatures: {
    preferred: 'ed25519',
    alternatives: ['secp256k1', 'p256'],
    notes: 'Ed25519 for speed, secp256k1 for Bitcoin/ETH compatibility, P-256 for WebAuthn',
  },
  pairing: {
    preferred: 'bls12-381',
    alternatives: ['bn254'],
    notes: 'BLS12-381 has better security margin; BN254 only if Ethereum precompiles needed',
  },
  zkCircuits: {
    preferred: 'bn254',
    alternatives: ['bls12-381', 'pallas'],
    notes: 'BN254 for Ethereum on-chain verification; BLS12-381/Pallas for higher security',
  },
  recursiveProofs: {
    preferred: 'pasta',
    alternatives: ['bn254+grumpkin'],
    notes: 'Pasta (Pallas/Vesta) for Halo2; BN254+Grumpkin for Noir',
  },
};

// Embedding degrees for pairing curves (affects security)
export const embeddingDegrees: Record<string, number> = {
  bn254: 12,
  bls12_381: 12,
};

export const security = {
  levels: securityLevels,
  recommendations,
  embeddingDegrees,
};
