/**
 * # Applications of Pairings in Cryptography
 *
 * Pairings enable cryptographic constructions that were previously thought
 * impossible. This tutorial explores the major applications beyond BLS signatures.
 *
 * ## Overview of Applications
 *
 * 1. Identity-Based Encryption (IBE)
 * 2. Attribute-Based Encryption (ABE)
 * 3. Short Signatures (beyond BLS)
 * 4. Zero-Knowledge Proofs (SNARKs)
 * 5. Verifiable Random Functions (VRF)
 * 6. Key Exchange Protocols
 *
 * ## The Unifying Theme
 *
 * Pairings let us verify relationships between discrete logs:
 *   - Given g^a and g^b, check if some value equals g^(ab)
 *
 * This simple capability enables remarkable constructions.
 *
 * @module tutorial/part6-advanced/17-pairings/pairing-applications
 */

import { GF, EllipticCurve } from 'sagemath-ts';

console.log("=== Pairing Applications in Cryptography ===\n");

// ============================================================================
// Section 1: Identity-Based Encryption (IBE)
// ============================================================================

console.log("--- Section 1: Identity-Based Encryption ---\n");

/**
 * IDENTITY-BASED ENCRYPTION (Boneh-Franklin 2001)
 *
 * Traditional PKI Problem:
 * To send encrypted email to alice@example.com, you need:
 * 1. Look up Alice's public key certificate
 * 2. Verify the certificate chain
 * 3. Check for revocation
 * 4. Only then can you encrypt!
 *
 * IBE Solution:
 * Alice's public key IS her identity (email address).
 * Anyone can encrypt to "alice@example.com" directly!
 *
 * SETUP:
 * - Key Generation Center (KGC) with master key msk
 * - Public parameters (pairing groups, hash function, mpk = msk * g2)
 *
 * KEY EXTRACTION:
 * - Alice authenticates to KGC
 * - KGC computes: sk_Alice = msk * H(alice@example.com)
 * - Alice receives her private key over secure channel
 *
 * ENCRYPTION (to identity ID):
 * - Sender computes: Q_ID = H(ID) in G1
 * - Choose random r
 * - Ciphertext: (r*g, m XOR H2(e(Q_ID, mpk)^r))
 *
 * DECRYPTION:
 * - Compute: e(C1, sk_ID) = e(r*g, msk*H(ID)) = e(g, H(ID))^(r*msk)
 *                        = e(H(ID), mpk)^r
 * - Recover m by XORing with H2(...)
 *
 * ADVANTAGES:
 * - No certificate lookup needed
 * - Natural key escrow (KGC can decrypt)
 * - Ideal for enterprise/closed systems
 *
 * LIMITATIONS:
 * - KGC has master key (key escrow by design)
 * - Requires secure channel for key extraction
 */

console.log("Identity-Based Encryption (Boneh-Franklin):\n");

console.log("Traditional PKI Problem:");
console.log("  To encrypt to Alice, must find and verify her certificate\n");

console.log("IBE Solution:");
console.log("  Public key = Identity (email, name, etc.)");
console.log("  Encrypt directly to 'alice@example.com'!\n");

console.log("Key Extraction:");
console.log("  sk_ID = msk * H(identity)");
console.log("  Only works because pairings enable the decryption math\n");

console.log("Use cases:");
console.log("  - Enterprise email encryption");
console.log("  - Certificate-less secure messaging");
console.log("  - Time-locked encryption (future timestamps as IDs)\n");

// ============================================================================
// Section 2: Attribute-Based Encryption (ABE)
// ============================================================================

console.log("--- Section 2: Attribute-Based Encryption ---\n");

/**
 * ATTRIBUTE-BASED ENCRYPTION
 *
 * Generalization of IBE where keys and ciphertexts are associated with
 * attributes and policies.
 *
 * KEY-POLICY ABE (KP-ABE):
 * - Ciphertext tagged with attributes: {doctor, cardiology, hospital-A}
 * - Secret key associated with policy: (doctor AND hospital-A) OR admin
 * - Decryption succeeds if attributes satisfy policy
 *
 * CIPHERTEXT-POLICY ABE (CP-ABE):
 * - Ciphertext has policy: (clearance=TOP_SECRET AND department=RESEARCH)
 * - Secret key has attributes: {clearance=TOP_SECRET, dept=RESEARCH, ...}
 * - Decryption succeeds if attributes satisfy ciphertext's policy
 *
 * CONSTRUCTION (simplified):
 * - Master key generates attribute keys
 * - Encryption embeds policy using Lagrange interpolation
 * - Decryption requires combining correct attribute keys
 * - Pairing enables checking attribute combinations
 *
 * EXAMPLE USE CASE (Healthcare):
 *
 * Policy: "(role=doctor AND dept=oncology) OR role=patient_self"
 *
 * Dr. Smith's key: {role=doctor, dept=oncology, employee_id=12345}
 *   -> Can decrypt medical records for oncology
 *
 * Patient's key: {role=patient_self, patient_id=67890}
 *   -> Can decrypt their own records
 *
 * Nurse's key: {role=nurse, dept=oncology}
 *   -> CANNOT decrypt (nurse != doctor)
 *
 * ADVANCED ABE:
 * - Multi-authority ABE: No single trusted party
 * - Decentralized ABE: Collusion-resistant even across authorities
 * - Revocable ABE: Can revoke specific user's access
 */

console.log("Attribute-Based Encryption (ABE):\n");

console.log("Key-Policy ABE:");
console.log("  Ciphertext: tagged with attributes");
console.log("  Secret key: embeds access policy");
console.log("  Example: Encrypt with {secret, project-X}");
console.log("           Decrypt if key has (secret AND project-X) OR admin\n");

console.log("Ciphertext-Policy ABE:");
console.log("  Ciphertext: contains access policy");
console.log("  Secret key: has user's attributes");
console.log("  Example: Encrypt with policy (clearance>=SECRET)");
console.log("           Decrypt if user has clearance attribute\n");

console.log("Healthcare Example:");
console.log("  Policy: (doctor AND oncology) OR patient_self");
console.log("  Dr. Smith (oncology):    CAN decrypt");
console.log("  Patient (own records):   CAN decrypt");
console.log("  Nurse (not doctor):      CANNOT decrypt\n");

// ============================================================================
// Section 3: Zero-Knowledge Proofs (SNARKs)
// ============================================================================

console.log("--- Section 3: SNARKs (Zero-Knowledge Proofs) ---\n");

/**
 * SNARKs: SUCCINCT NON-INTERACTIVE ARGUMENTS OF KNOWLEDGE
 *
 * SNARKs are perhaps the most impactful application of pairings today.
 * They enable proving arbitrary computations with:
 * - Succinct proofs: ~200-300 bytes regardless of computation size
 * - Fast verification: milliseconds regardless of computation
 * - Zero-knowledge: proof reveals nothing about inputs
 *
 * WHY PAIRINGS?
 *
 * SNARKs work by encoding computations as polynomial equations.
 * The prover must show they know polynomials satisfying constraints
 * WITHOUT revealing the polynomials.
 *
 * Pairings enable polynomial commitment schemes:
 * - Commit to polynomial f(x) as [f(s)]_1 for secret s
 * - Prove f(z) = y without revealing f
 * - Uses: e([f(s)]_1, [1]_2) checks
 *
 * GROTH16 (most common SNARK):
 * - Setup: Trusted ceremony generates proving/verification keys
 * - Prove: ~seconds for complex computations
 * - Verify: ~milliseconds, 3 pairings
 * - Proof size: 192 bytes (3 G1 elements)
 *
 * PLONK:
 * - Universal trusted setup (one ceremony for all circuits)
 * - More flexible constraint system
 * - Slightly larger proofs than Groth16
 *
 * APPLICATIONS:
 * - Zcash: Private transactions (prove valid without revealing amounts)
 * - Ethereum: Rollups (prove thousands of txs in one proof)
 * - Filecoin: Storage proofs (prove you're storing data)
 * - Identity: Prove age > 18 without revealing birthday
 */

console.log("SNARKs - Zero-Knowledge Proofs:\n");

console.log("Properties:");
console.log("  - Succinct: ~200 bytes regardless of computation");
console.log("  - Non-interactive: no back-and-forth");
console.log("  - Zero-knowledge: reveals nothing about inputs\n");

console.log("Why Pairings are Essential:");
console.log("  - Polynomial commitment: commit to f(x) as [f(s)]_1");
console.log("  - Pairing checks: e([f(s)]_1, [g]_2) enables verification");
console.log("  - KZG commitments use pairing for opening proofs\n");

console.log("Groth16 (most efficient):");
console.log("  - Proof: 192 bytes (3 G1 points)");
console.log("  - Verify: 3 pairings, ~1ms");
console.log("  - Trusted setup per circuit\n");

console.log("Applications:");
console.log("  - Zcash: Private transactions");
console.log("  - Rollups: Scale blockchains 100x+");
console.log("  - Privacy: Prove attributes without revealing data\n");

// ============================================================================
// Section 4: Verifiable Random Functions (VRF)
// ============================================================================

console.log("--- Section 4: Verifiable Random Functions ---\n");

/**
 * VERIFIABLE RANDOM FUNCTIONS (VRF)
 *
 * A VRF is like a keyed hash where you can PROVE the output is correct.
 *
 * PROPERTIES:
 * - Deterministic: VRF(sk, x) always gives same output
 * - Unpredictable: Without sk, output looks random
 * - Verifiable: Can prove output is correct using pk
 *
 * PAIRING-BASED VRF:
 *
 * KeyGen:
 *   sk <- random, pk = sk * g
 *
 * Eval(sk, x):
 *   h = H(x)               // Hash to curve
 *   gamma = sk * h          // VRF "proof point"
 *   output = H'(gamma)      // Final random output
 *   return (output, gamma)
 *
 * Verify(pk, x, output, gamma):
 *   h = H(x)
 *   Check: e(gamma, g) == e(h, pk)  // Pairing verification!
 *   Check: output == H'(gamma)
 *
 * WHY IT WORKS:
 *   e(sk * h, g) = e(h, g)^sk = e(h, sk * g) = e(h, pk)
 *
 * The pairing proves gamma = sk * H(x) without revealing sk!
 *
 * APPLICATIONS:
 * - Blockchain leader election (Cardano, Algorand)
 * - Random beacons (unbiasable randomness)
 * - Lottery systems (prove fairness)
 * - Private information retrieval
 */

console.log("Verifiable Random Functions (VRF):\n");

console.log("What it provides:");
console.log("  - Deterministic random-looking output");
console.log("  - Proof that output is correct");
console.log("  - Without revealing secret key\n");

console.log("Construction:");
console.log("  gamma = sk * H(input)");
console.log("  output = H'(gamma)");
console.log("  Verify: e(gamma, g) ?= e(H(input), pk)\n");

console.log("Applications:");
console.log("  - Proof-of-Stake leader election");
console.log("  - Unbiasable random beacons");
console.log("  - Provably fair lotteries\n");

// ============================================================================
// Section 5: Tripartite Key Exchange
// ============================================================================

console.log("--- Section 5: Key Exchange Protocols ---\n");

/**
 * TRIPARTITE KEY EXCHANGE (Joux 2000)
 *
 * Classical DH requires two rounds for 3 parties.
 * Pairings enable ONE-ROUND tripartite key exchange!
 *
 * PROTOCOL:
 * - Alice, Bob, Carol each pick secrets a, b, c
 * - Alice broadcasts aP
 * - Bob broadcasts bP
 * - Carol broadcasts cP
 *
 * Shared key computation:
 * - Alice: K = e(bP, cP)^a
 * - Bob:   K = e(aP, cP)^b
 * - Carol: K = e(aP, bP)^c
 *
 * All compute the same K = e(P, P)^{abc}!
 *
 * PROOF (bilinearity):
 *   e(bP, cP)^a = e(P, P)^{bc * a} = e(P, P)^{abc}
 *
 * SECURITY:
 * - Based on Bilinear Diffie-Hellman (BDH) assumption
 * - Passive adversary cannot compute e(P, P)^{abc}
 *   from (P, aP, bP, cP)
 *
 * Note: This basic protocol lacks authentication.
 * Practical protocols add signatures or MACs.
 *
 * IDENTITY-BASED KEY EXCHANGE:
 * Combine IBE with key exchange for certificate-free authenticated
 * key exchange. Users derive keys from identities!
 */

console.log("Tripartite Key Exchange (Joux):\n");

console.log("Setup: Alice(a), Bob(b), Carol(c)");
console.log("  Each broadcasts their public point\n");

console.log("Key Derivation (one round!):");
console.log("  Alice: K = e(bP, cP)^a");
console.log("  Bob:   K = e(aP, cP)^b");
console.log("  Carol: K = e(aP, bP)^c");
console.log("  All get: K = e(P, P)^{abc}\n");

console.log("Comparison:");
console.log("  Traditional 3-party DH: 2 rounds");
console.log("  Pairing-based: 1 round!");
console.log("  Security: BDH assumption\n");

// ============================================================================
// Section 6: More Applications
// ============================================================================

console.log("--- Section 6: Additional Applications ---\n");

/**
 * MORE PAIRING APPLICATIONS
 *
 * SHORT SIGNATURES (Waters, etc.):
 * - Even shorter than BLS in some settings
 * - Structure-preserving signatures
 * - Useful for anonymous credentials
 *
 * GROUP SIGNATURES:
 * - User signs on behalf of group
 * - Verifier knows "someone in group signed"
 * - Group manager can open to reveal signer
 * - Privacy + accountability
 *
 * BLIND SIGNATURES:
 * - Signer signs without seeing message
 * - Used for anonymous voting, e-cash
 * - Pairing-based schemes very efficient
 *
 * FUNCTIONAL ENCRYPTION:
 * - Decrypt only f(m) from encrypted m
 * - Example: Compute on encrypted data
 * - Inner-product functional encryption from pairings
 *
 * PROXY RE-ENCRYPTION:
 * - Transform ciphertext from pk_A to pk_B
 * - Without decrypting or learning message
 * - Useful for secure data sharing
 *
 * BROADCAST ENCRYPTION:
 * - Encrypt to subset of users efficiently
 * - Revoke users without re-keying others
 * - Pairing-based schemes achieve optimal parameters
 */

console.log("Additional Applications:\n");

console.log("Group Signatures:");
console.log("  Sign on behalf of group");
console.log("  Verify: 'someone in group signed'");
console.log("  Open: manager can identify signer\n");

console.log("Functional Encryption:");
console.log("  Decrypt f(m) without learning m");
console.log("  Compute on encrypted data\n");

console.log("Proxy Re-Encryption:");
console.log("  Transform Enc(pk_A, m) to Enc(pk_B, m)");
console.log("  Without decrypting!\n");

console.log("Broadcast Encryption:");
console.log("  Encrypt to subset of users");
console.log("  Efficiently revoke access\n");

// ============================================================================
// Section 7: Choosing the Right Application
// ============================================================================

console.log("--- Section 7: When to Use Pairings ---\n");

/**
 * WHEN TO USE PAIRING-BASED CRYPTOGRAPHY
 *
 * GOOD USE CASES:
 *
 * 1. Signature Aggregation:
 *    - Many signers, limited bandwidth
 *    - Example: Blockchain consensus, IoT sensor networks
 *
 * 2. Identity-Based Systems:
 *    - Controlled environments with key authority
 *    - Example: Enterprise encryption, smart cards
 *
 * 3. Advanced Access Control:
 *    - Complex policies, attribute-based
 *    - Example: Healthcare records, classified documents
 *
 * 4. Zero-Knowledge Proofs:
 *    - Proving computations without revealing inputs
 *    - Example: Privacy coins, scaling solutions
 *
 * 5. Threshold Cryptography:
 *    - Distributed trust, no single point of failure
 *    - Example: Key custody, random beacons
 *
 * LESS SUITABLE:
 *
 * 1. Simple Encryption:
 *    - ECIES or hybrid encryption is simpler
 *
 * 2. Resource-Constrained Devices:
 *    - Pairing computation is heavy
 *    - Consider ECDSA for IoT
 *
 * 3. Post-Quantum Security:
 *    - Pairings are NOT quantum-resistant
 *    - Use lattice-based alternatives for PQ
 *
 * PERFORMANCE CONSIDERATIONS:
 * - Pairing: ~1-2ms (BLS12-381)
 * - EC scalar mult: ~0.1-0.2ms
 * - RSA-2048: ~0.5ms sign, ~0.02ms verify
 * - Post-quantum (Dilithium): ~0.1ms
 */

console.log("Good Use Cases for Pairings:\n");
console.log("  - Signature aggregation (BLS)");
console.log("  - Identity-based encryption");
console.log("  - Attribute-based encryption");
console.log("  - Zero-knowledge proofs (SNARKs)");
console.log("  - Threshold cryptography\n");

console.log("Consider Alternatives When:\n");
console.log("  - Simple encryption needed (use ECIES)");
console.log("  - Very constrained devices (use ECDSA)");
console.log("  - Post-quantum required (use lattices)\n");

console.log("Performance Comparison (BLS12-381):");
console.log("  Pairing:      ~1-2ms");
console.log("  EC mult:      ~0.1-0.2ms");
console.log("  BLS sign:     ~0.5ms");
console.log("  BLS verify:   ~1.5ms (includes pairing)");
console.log("  Aggregation:  Amortizes verification cost\n");

// ============================================================================
// Summary
// ============================================================================

console.log("=== Summary ===\n");

console.log("Pairings enable cryptographic 'magic':\n");

console.log("Identity-Based Encryption:");
console.log("  Public key = email address, no certificates!\n");

console.log("Attribute-Based Encryption:");
console.log("  Access control embedded in encryption\n");

console.log("SNARKs:");
console.log("  Prove any computation in ~200 bytes\n");

console.log("Verifiable Random Functions:");
console.log("  Prove randomness is correctly generated\n");

console.log("Key Exchange:");
console.log("  3-party in one round (impossible without pairings)\n");

console.log("The common thread: pairings let us verify");
console.log("relationships between discrete logarithms.\n");

console.log("This concludes the pairings chapter.");
console.log("Next: Introduction to Lattice Cryptography");
