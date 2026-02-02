/**
 * # Lattice-Based Cryptography: The Post-Quantum Future
 *
 * Lattice-based cryptography is the leading approach to post-quantum security.
 * NIST has standardized lattice-based schemes for encryption (Kyber) and
 * signatures (Dilithium), which will replace RSA and ECC in coming years.
 *
 * ## Why Lattices for Post-Quantum?
 *
 * 1. NO KNOWN QUANTUM ATTACKS: Unlike factoring/DLP, no quantum speedup
 * 2. EFFICIENT: Often faster than RSA, competitive with ECC
 * 3. VERSATILE: Enables FHE, IBE, ABE, ZK proofs
 * 4. WELL-STUDIED: Decades of cryptanalysis, strong reductions
 *
 * This tutorial covers the key concepts: LWE, Ring-LWE, and the schemes
 * built from them.
 *
 * @module tutorial/part6-advanced/18-lattices-intro/lattice-crypto
 */

import { LWE, Regev, LindnerPeikert, UniformSampler, samples } from 'sagemath-ts';

console.log("=== Lattice-Based Cryptography ===\n");

// ============================================================================
// Section 1: Learning With Errors (LWE)
// ============================================================================

console.log("--- Section 1: Learning With Errors (LWE) ---\n");

/**
 * LEARNING WITH ERRORS (LWE)
 *
 * The LWE problem is the foundation of most lattice cryptography.
 *
 * SETUP:
 * - Dimension n (security parameter)
 * - Modulus q (typically polynomial in n)
 * - Error distribution chi (typically Gaussian with small std dev)
 * - Secret vector s in Z_q^n
 *
 * LWE DISTRIBUTION:
 * Generate samples (a, b) where:
 *   a <- uniform random in Z_q^n
 *   e <- error from chi
 *   b = <a, s> + e mod q
 *
 * SEARCH LWE:
 * Given many samples (a_i, b_i), find s
 *
 * DECISION LWE:
 * Distinguish LWE samples from uniform random pairs
 *
 * WHY IT'S HARD:
 *
 * Without noise: Solving linear systems is easy (Gaussian elimination)
 * With noise:    The error "hides" the linear relationship
 *
 * SECURITY REDUCTION (Regev 2005):
 * If there's an efficient algorithm for LWE,
 * there's an efficient quantum algorithm for worst-case lattice problems
 * (GapSVP, SIVP with gamma = O(n/alpha))
 *
 * This is a WORST-CASE TO AVERAGE-CASE reduction!
 * Breaking random LWE instances implies solving any lattice instance.
 */

console.log("LWE (Learning With Errors):\n");
console.log("  Setup: dimension n, modulus q, error distribution chi, secret s\n");
console.log("  Sample: a <- random, e <- chi, b = <a,s> + e mod q\n");
console.log("  Problem: Given many (a, b) pairs, find s\n");

console.log("Why it's hard:");
console.log("  Without noise: Easy (linear algebra)");
console.log("  With noise:    Noise hides linear relationship\n");

console.log("Security:");
console.log("  Worst-case to average-case reduction");
console.log("  Breaking LWE => Solving worst-case SVP/SIVP\n");

// Demo with sagemath-ts LWE
console.log("LWE Demo:\n");

const lweN = 20n;  // Small for demo (real crypto uses ~512-1024)
const lwe = new Regev(lweN);

console.log(`  Parameters: n = ${lweN}, q = ${lwe.K.modulus}`);
console.log(`  Secret: (hidden, but has ${lweN} components)`);

// Generate some samples
const [a1, b1] = lwe.call();
const [a2, b2] = lwe.call();

console.log("\n  Sample 1:");
console.log(`    a = (${a1.entries.slice(0, 5).map(x => x.value).join(', ')}, ...)`);
console.log(`    b = ${b1.value}`);

console.log("\n  Sample 2:");
console.log(`    a = (${a2.entries.slice(0, 5).map(x => x.value).join(', ')}, ...)`);
console.log(`    b = ${b2.value}\n`);

// ============================================================================
// Section 2: Ring-LWE and Module-LWE
// ============================================================================

console.log("--- Section 2: Ring-LWE and Module-LWE ---\n");

/**
 * STRUCTURED LWE VARIANTS
 *
 * Plain LWE has large keys: secret s has n independent components.
 * For n = 1024, keys are ~1024 * log(q) bits.
 *
 * RING-LWE:
 * Work in polynomial ring R_q = Z_q[x] / (x^n + 1) for n = power of 2.
 *
 * - Secret s is now ONE polynomial in R_q
 * - Sample a is also ONE polynomial
 * - Product a*s is polynomial multiplication (fast via NTT)
 * - Keys are n times smaller!
 *
 * Sample: a <- R_q, e <- small polynomial, b = a*s + e in R_q
 *
 * SECURITY:
 * Reduction to worst-case ideal lattice problems.
 * Slightly stronger assumption than plain LWE.
 *
 * MODULE-LWE:
 * Interpolation between LWE and Ring-LWE.
 * Work with k x k matrices of ring elements.
 *
 * - k = 1: Ring-LWE
 * - k = n: Plain LWE (roughly)
 *
 * Module-LWE with k = 2, 3, or 4 gives good balance.
 * KYBER uses Module-LWE with k in {2, 3, 4}.
 *
 * ADVANTAGES OF RING/MODULE:
 * - Smaller keys and ciphertexts
 * - Fast operations (NTT for multiplication)
 * - Still quantum-resistant
 *
 * CONCERNS:
 * - Algebraic structure might help attacks
 * - No known concrete attacks, but less studied than plain LWE
 */

console.log("Ring-LWE:");
console.log("  Work in R_q = Z_q[x] / (x^n + 1)");
console.log("  Secret s is a polynomial (n coefficients)");
console.log("  Multiplication is polynomial mult (NTT: O(n log n))");
console.log("  Keys ~n times smaller than plain LWE\n");

console.log("Module-LWE (used in Kyber):");
console.log("  k x k matrices of ring elements");
console.log("  k = 2: Kyber-512  (~128-bit security)");
console.log("  k = 3: Kyber-768  (~192-bit security)");
console.log("  k = 4: Kyber-1024 (~256-bit security)\n");

console.log("Efficiency comparison (key sizes):");
console.log("  Plain LWE (n=1024): ~12 KB public key");
console.log("  Module-LWE (k=3, n=256): ~1.1 KB public key");
console.log("  Ring-LWE (n=1024): ~1 KB public key\n");

// ============================================================================
// Section 3: Kyber Key Encapsulation
// ============================================================================

console.log("--- Section 3: Kyber (NIST Standard) ---\n");

/**
 * KYBER - ML-KEM (Module-Lattice Key Encapsulation)
 *
 * Kyber is the NIST-standardized post-quantum key encapsulation mechanism.
 * It's based on Module-LWE.
 *
 * KEY GENERATION:
 * 1. Generate random secret s (vector of small polynomials)
 * 2. Generate random A (matrix of polynomials, from seed)
 * 3. Compute b = A*s + e for error e
 * 4. Public key: (A, b) [actually just seed + b]
 * 5. Secret key: s
 *
 * ENCAPSULATION (encrypt random key):
 * 1. Sample random r (used for encryption randomness)
 * 2. Compute u = A^T * r + e1
 * 3. Compute v = b^T * r + e2 + encode(key)
 * 4. Ciphertext: (u, v)
 *
 * DECAPSULATION (decrypt to get key):
 * 1. Compute v - s^T * u = (b^T - s^T * A^T) * r + noise + encode(key)
 *                       = (s^T * A + e^T - s^T * A^T) * r + noise + encode(key)
 *                       = e^T * r + noise + encode(key)
 * 2. The error is small, so round to decode the key
 *
 * PARAMETERS (Kyber-768):
 * - n = 256 (ring dimension)
 * - k = 3 (module rank)
 * - q = 3329 (modulus)
 * - Error distribution: binomial with eta = 2
 *
 * PERFORMANCE:
 * - KeyGen: ~30 microseconds
 * - Encaps: ~40 microseconds
 * - Decaps: ~40 microseconds
 * - Public key: 1,184 bytes
 * - Ciphertext: 1,088 bytes
 *
 * Much faster than RSA, comparable to ECC!
 */

console.log("Kyber (ML-KEM) - NIST Standard:\n");

console.log("Key Generation:");
console.log("  1. Sample secret s (vector of small polynomials)");
console.log("  2. Compute b = A*s + e");
console.log("  3. Public key: (seed for A, b)");
console.log("  4. Secret key: s\n");

console.log("Encapsulation (encrypt random key K):");
console.log("  1. Sample random r");
console.log("  2. u = A^T * r + e1");
console.log("  3. v = b^T * r + e2 + encode(K)");
console.log("  4. Ciphertext: (u, v)\n");

console.log("Decapsulation:");
console.log("  1. Compute v - s^T * u");
console.log("  2. Round to decode K\n");

console.log("Kyber-768 Parameters:");
console.log("  n = 256, k = 3, q = 3329");
console.log("  Public key: 1,184 bytes");
console.log("  Ciphertext: 1,088 bytes");
console.log("  ~128-bit post-quantum security\n");

console.log("Performance (compared to classical):");
console.log("  Kyber:  ~40us encaps, 1KB keys");
console.log("  RSA-2048: ~100us sign, 256B keys");
console.log("  ECC-256: ~50us ECDH, 32B keys");
console.log("  Kyber is competitive with classical!\n");

// ============================================================================
// Section 4: Dilithium Signatures
// ============================================================================

console.log("--- Section 4: Dilithium Signatures ---\n");

/**
 * DILITHIUM - ML-DSA (Module-Lattice Digital Signature)
 *
 * Dilithium is the NIST-standardized post-quantum signature scheme.
 * Based on the "Fiat-Shamir with aborts" technique.
 *
 * BASED ON:
 * Short Integer Solution (SIS) problem:
 *   Given random A, find short z with A*z = 0 mod q
 *
 * KEY GENERATION:
 * 1. Generate random A (from seed)
 * 2. Sample short secrets s1, s2
 * 3. Compute t = A*s1 + s2
 * 4. Public key: (A, t) [seed + t]
 * 5. Secret key: (s1, s2)
 *
 * SIGNING:
 * 1. Sample random y (masking vector)
 * 2. Compute w = A*y
 * 3. Compute challenge c = H(w, message)
 * 4. Compute z = y + c*s1
 * 5. If z is too large, REJECT and try again (abort)
 * 6. Signature: (z, h) where h encodes some hint
 *
 * VERIFICATION:
 * 1. Compute challenge c = H(A*z - c*t, message)
 * 2. Check that z is small enough
 *
 * WHY ABORT?
 * Without rejection sampling, z would leak information about s1.
 * By rejecting "bad" z values, the distribution of z becomes
 * independent of s1.
 *
 * PARAMETERS (Dilithium3):
 * - n = 256, k = 6, l = 5
 * - q = 8380417
 * - Signature: 3,293 bytes
 * - Public key: 1,952 bytes
 * - ~128-bit security
 */

console.log("Dilithium (ML-DSA) - NIST Standard:\n");

console.log("Based on SIS problem:");
console.log("  Find short z with A*z = 0 mod q\n");

console.log("Key Generation:");
console.log("  1. Sample short s1, s2");
console.log("  2. Compute t = A*s1 + s2");
console.log("  3. Public key: (A, t)\n");

console.log("Signing (with rejection sampling):");
console.log("  1. Sample random y, compute w = A*y");
console.log("  2. Challenge c = H(w, message)");
console.log("  3. z = y + c*s1");
console.log("  4. If ||z|| too large, RESTART");
console.log("  5. Signature: (z, hint)\n");

console.log("Verification:");
console.log("  1. Recompute challenge from A*z - c*t");
console.log("  2. Verify z is small\n");

console.log("Dilithium3 Parameters:");
console.log("  Signature: 3,293 bytes");
console.log("  Public key: 1,952 bytes");
console.log("  Private key: 4,016 bytes");
console.log("  ~128-bit post-quantum security\n");

// ============================================================================
// Section 5: NTRU
// ============================================================================

console.log("--- Section 5: NTRU (Alternative Lattice Scheme) ---\n");

/**
 * NTRU - OLDER LATTICE ENCRYPTION
 *
 * NTRU predates LWE (published 1998) and has a different structure.
 * It's also being standardized as an alternative to Kyber.
 *
 * RING:
 * Work in R = Z[x] / (x^n - 1) (note: NOT x^n + 1)
 * Elements are polynomials of degree < n
 *
 * KEY GENERATION:
 * 1. Sample small polynomials f, g
 * 2. Require f invertible mod q and mod p
 * 3. Public key h = g * f^{-1} mod q
 * 4. Secret key: f (and derived values)
 *
 * ENCRYPTION:
 * 1. Sample small random r
 * 2. Ciphertext c = p*r*h + m mod q
 *
 * DECRYPTION:
 * 1. Compute a = f * c mod q
 *            = f * (p*r*h + m) mod q
 *            = f * (p*r*g*f^{-1} + m) mod q
 *            = p*r*g + f*m mod q
 * 2. Reduce a mod p to get f*m mod p
 * 3. Multiply by f^{-1} mod p to get m
 *
 * This works because p*r*g is "small" and reduces to 0 mod p.
 *
 * COMPARISON WITH KYBER:
 * - NTRU: older, more studied algebraic structure
 * - Kyber: based on (Module-)LWE, cleaner security proofs
 * - Both are being standardized by NIST
 *
 * NTRU is also efficient and has very compact ciphertexts.
 */

console.log("NTRU (1998):\n");

console.log("Ring: Z[x] / (x^n - 1)");
console.log("  Different from LWE-based (which use x^n + 1)\n");

console.log("Key Generation:");
console.log("  1. Sample small f, g");
console.log("  2. Public key h = g * f^{-1} mod q");
console.log("  3. Secret key: f\n");

console.log("Encryption of m:");
console.log("  c = p*r*h + m mod q\n");

console.log("Decryption:");
console.log("  1. Compute f*c = p*r*g + f*m mod q");
console.log("  2. Reduce mod p (small error vanishes)");
console.log("  3. Multiply by f^{-1} mod p\n");

console.log("NTRU vs Kyber:");
console.log("  NTRU:  Older, well-studied, compact");
console.log("  Kyber: Better security proofs, modular");
console.log("  Both are NIST PQC candidates\n");

// ============================================================================
// Section 6: Fully Homomorphic Encryption
// ============================================================================

console.log("--- Section 6: FHE (Fully Homomorphic Encryption) ---\n");

/**
 * FULLY HOMOMORPHIC ENCRYPTION (FHE)
 *
 * FHE allows computation on encrypted data without decryption.
 * The holy grail of cryptography, first achieved in 2009 using lattices!
 *
 * WHAT FHE ENABLES:
 * - Cloud computing on private data
 * - Privacy-preserving machine learning
 * - Secure multi-party computation
 * - Private database queries
 *
 * LWE-BASED FHE:
 *
 * Encryption: ct = (a, b = <a,s> + e + m*q/2) for bit m
 *
 * ADDITION: ct1 + ct2 = (a1+a2, b1+b2)
 *   Decrypts to m1 + m2 (errors add up)
 *
 * MULTIPLICATION: Much more complex!
 *   - Requires "relinearization" to keep ciphertext size bounded
 *   - Noise grows multiplicatively
 *
 * NOISE MANAGEMENT (The Key Challenge):
 * - Each operation increases noise
 * - Eventually noise exceeds threshold and decryption fails
 * - BOOTSTRAPPING: Homomorphically evaluate decryption to reduce noise
 * - Bootstrapping enables unlimited computation (but slow!)
 *
 * SCHEMES:
 * - BGV/BFV: Integer arithmetic, good for ML
 * - CKKS: Approximate arithmetic, good for floats
 * - TFHE: Boolean circuits, fastest bootstrapping
 *
 * PERFORMANCE:
 * - Still 1000-10000x slower than plaintext
 * - Active research area
 * - Practical for some applications today
 */

console.log("Fully Homomorphic Encryption (FHE):\n");

console.log("The dream: Compute on encrypted data!");
console.log("  Enc(a) + Enc(b) = Enc(a + b)");
console.log("  Enc(a) * Enc(b) = Enc(a * b)\n");

console.log("Why lattices enable FHE:");
console.log("  LWE ciphertexts support natural addition");
console.log("  Multiplication is harder (noise grows)");
console.log("  Bootstrapping resets noise (key breakthrough)\n");

console.log("FHE Schemes:");
console.log("  BGV/BFV: Integer arithmetic");
console.log("  CKKS:    Approximate (for ML on floats)");
console.log("  TFHE:    Boolean gates (fast bootstrapping)\n");

console.log("Current state:");
console.log("  ~1000-10000x overhead vs plaintext");
console.log("  Practical for some applications");
console.log("  Active research to improve performance\n");

// ============================================================================
// Section 7: The Post-Quantum Transition
// ============================================================================

console.log("--- Section 7: The Post-Quantum Transition ---\n");

/**
 * MIGRATION TO POST-QUANTUM CRYPTOGRAPHY
 *
 * NIST PQC STANDARDS (2024):
 *
 * Key Encapsulation:
 * - ML-KEM (Kyber): Primary standard
 * - Alternate candidates in evaluation
 *
 * Digital Signatures:
 * - ML-DSA (Dilithium): Primary lattice-based
 * - SLH-DSA (SPHINCS+): Hash-based backup
 * - More signatures being evaluated
 *
 * TIMELINE:
 * 2024: Final standards published
 * 2025-2030: Transition period
 * 2030+: Deprecation of RSA/ECC for sensitive data
 *
 * "HARVEST NOW, DECRYPT LATER" THREAT:
 * Adversaries may store encrypted data today,
 * wait for quantum computers, then decrypt.
 * This motivates EARLY migration for sensitive data.
 *
 * HYBRID APPROACH:
 * Combine classical + PQC algorithms during transition:
 * - TLS: Use both ECDH and Kyber
 * - Signatures: Use both ECDSA and Dilithium
 * If either breaks, the other provides security.
 *
 * IMPLEMENTATION CHALLENGES:
 * - Larger keys and ciphertexts (network overhead)
 * - Different side-channel characteristics
 * - New implementation bugs to discover
 * - Hardware acceleration not yet widespread
 */

console.log("NIST PQC Standards (2024):\n");
console.log("  ML-KEM (Kyber):     Key encapsulation");
console.log("  ML-DSA (Dilithium): Digital signatures");
console.log("  SLH-DSA (SPHINCS+): Hash-based signatures\n");

console.log("Timeline:");
console.log("  2024:      Standards published");
console.log("  2025-2030: Transition period");
console.log("  2030+:     RSA/ECC deprecated for sensitive data\n");

console.log("Hybrid Approach (recommended for transition):");
console.log("  Combine: ECDH + Kyber for key exchange");
console.log("  Combine: ECDSA + Dilithium for signatures");
console.log("  Security if EITHER is secure\n");

console.log("'Harvest now, decrypt later' threat:");
console.log("  Adversaries storing data for future quantum attack");
console.log("  Reason to migrate NOW for long-term secrets\n");

// ============================================================================
// Summary
// ============================================================================

console.log("=== Summary ===\n");

console.log("Lattice-Based Cryptography provides post-quantum security:\n");

console.log("LWE Problem:");
console.log("  b = <a, s> + e is hard to solve");
console.log("  Basis of Kyber, FHE, and more\n");

console.log("Ring/Module-LWE:");
console.log("  Use polynomial rings for efficiency");
console.log("  Smaller keys, faster operations\n");

console.log("NIST Standards:");
console.log("  Kyber (ML-KEM): Post-quantum key exchange");
console.log("  Dilithium (ML-DSA): Post-quantum signatures\n");

console.log("Performance:");
console.log("  Competitive with classical crypto!");
console.log("  Larger keys but similar speed\n");

console.log("The post-quantum transition is happening NOW.");
console.log("Lattice cryptography is the leading solution.\n");

console.log("=== End of Part 6: Advanced Topics ===");
