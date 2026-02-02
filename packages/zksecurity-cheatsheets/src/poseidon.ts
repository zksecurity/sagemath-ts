export interface PoseidonParams {
  name: string;
  field: bigint; // Prime field modulus
  t: number; // State width (rate + capacity)
  rate: number; // Rate (number of field elements absorbed per permutation)
  capacity: number; // Capacity (security parameter, usually 1)
  fullRounds: number; // Number of full rounds (Rf)
  partialRounds: number; // Number of partial rounds (Rp)
  alpha: number; // S-box exponent (usually 5 or 7)
  description: string;
}

// BN254 Poseidon configurations (used by Circom, snarkjs)
export const poseidonBn254T3: PoseidonParams = {
  name: 'poseidon-bn254-t3',
  field: 21888242871839275222246405745257275088548364400416034343698204186575808495617n,
  t: 3,
  rate: 2,
  capacity: 1,
  fullRounds: 8,
  partialRounds: 57,
  alpha: 5,
  description: 'Circom/snarkjs default for 2-to-1 hash',
};

export const poseidonBn254T5: PoseidonParams = {
  name: 'poseidon-bn254-t5',
  field: 21888242871839275222246405745257275088548364400416034343698204186575808495617n,
  t: 5,
  rate: 4,
  capacity: 1,
  fullRounds: 8,
  partialRounds: 60,
  alpha: 5,
  description: 'Circom/snarkjs for 4-to-1 hash',
};

// BLS12-381 Poseidon (used by Filecoin)
export const poseidonBls12_381T3: PoseidonParams = {
  name: 'poseidon-bls12-381-t3',
  field: 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001n,
  t: 3,
  rate: 2,
  capacity: 1,
  fullRounds: 8,
  partialRounds: 55,
  alpha: 5,
  description: 'Filecoin/Neptune for 2-to-1 hash',
};

// Goldilocks Poseidon (used by Plonky2)
export const poseidonGoldilocks: PoseidonParams = {
  name: 'poseidon-goldilocks-t12',
  field: 2n ** 64n - 2n ** 32n + 1n,
  t: 12,
  rate: 8,
  capacity: 4,
  fullRounds: 8,
  partialRounds: 22,
  alpha: 7,
  description: 'Plonky2 sponge construction',
};

// Pallas Poseidon (used by Halo2)
export const poseidonPallas: PoseidonParams = {
  name: 'poseidon-pallas-t3',
  field: 0x40000000000000000000000000000000224698fc094cf91b992d30ed00000001n,
  t: 3,
  rate: 2,
  capacity: 1,
  fullRounds: 8,
  partialRounds: 56,
  alpha: 5,
  description: 'Halo2/Zcash Orchard',
};

export const poseidon = {
  bn254T3: poseidonBn254T3,
  bn254T5: poseidonBn254T5,
  bls12_381T3: poseidonBls12_381T3,
  goldilocks: poseidonGoldilocks,
  pallas: poseidonPallas,
};
