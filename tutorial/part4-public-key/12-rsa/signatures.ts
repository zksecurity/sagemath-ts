/**
 * # RSA Digital Signatures
 *
 * Digital signatures provide authentication, integrity, and non-repudiation.
 * RSA signatures work by "encrypting" with the private key.
 *
 * ## Signature Scheme
 *
 * - Sign:   s = H(m)^d mod n    (sign hash of message)
 * - Verify: H(m) = s^e mod n    (verify using public key)
 *
 * ## Why Signatures Work
 *
 * Only the holder of the private key d can compute s = h^d mod n.
 * Anyone with the public key (n, e) can verify s^e ≡ h (mod n).
 *
 * ## Security Considerations
 *
 * - Always sign a hash of the message, not the message itself
 * - Use proper padding schemes (PSS for RSA signatures)
 * - The hash function must be collision-resistant
 */

import { gcd, inverse_mod, power_mod, Mod } from 'sagemath-ts';

console.log('='.repeat(60));
console.log('RSA Digital Signatures Tutorial');
console.log('='.repeat(60));

// ============================================================================
// Setup: RSA Key Pair
// ============================================================================

console.log('\n=== RSA Key Setup ===\n');

const p = 61n;
const q = 53n;
const n = p * q;
const phi_n = (p - 1n) * (q - 1n);
const e = 17n;
const d = inverse_mod(e, phi_n);

console.log('Key Parameters:');
console.log(`  n = ${n} (public)`);
console.log(`  e = ${e} (public)`);
console.log(`  d = ${d} (private)`);

// ============================================================================
// Basic "Textbook" RSA Signature
// ============================================================================

console.log('\n=== Basic RSA Signature ===\n');

/**
 * Textbook RSA Signature (for education only, not secure!):
 *   Sign:   s = m^d mod n
 *   Verify: m' = s^e mod n, check if m' = m
 *
 * This is insecure because:
 * 1. Signatures are malleable
 * 2. Existential forgery is possible
 * 3. Long messages don't fit
 */

const message = 42n;
console.log(`Message to sign: m = ${message}`);

// Sign: compute s = m^d mod n
const signature = Mod(message, n).pow(d).value;
console.log(`\nSigning: s = m^d mod n`);
console.log(`  s = ${message}^${d} mod ${n}`);
console.log(`  s = ${signature}`);

// Verify: compute s^e mod n and check if it equals m
const verificationValue = Mod(signature, n).pow(e).value;
console.log(`\nVerification: m' = s^e mod n`);
console.log(`  m' = ${signature}^${e} mod ${n}`);
console.log(`  m' = ${verificationValue}`);
console.log(`  Signature valid: ${verificationValue === message}`);

// ============================================================================
// Why We Sign Hashes, Not Messages
// ============================================================================

console.log('\n=== Why We Sign Hashes ===\n');

/**
 * Problem 1: Message Size
 * RSA can only sign messages m < n. Real messages are much larger.
 *
 * Problem 2: Existential Forgery
 * Attacker can create valid signatures for "random" messages:
 * 1. Pick random s
 * 2. Compute m = s^e mod n
 * 3. (s, m) is a valid signature pair!
 *
 * Solution: Sign H(m) where H is a cryptographic hash function
 */

console.log('Problem: Existential Forgery Attack');
console.log('');

// Demonstrate the attack
const fakeSignature = 123n; // Attacker picks arbitrary signature
const forgedMessage = Mod(fakeSignature, n).pow(e).value;

console.log(`  1. Attacker chooses random s = ${fakeSignature}`);
console.log(`  2. Attacker computes m = s^e mod n = ${forgedMessage}`);
console.log(`  3. Verification: s^e mod n = ${forgedMessage} (equals m)`);
console.log(`  4. The pair (${fakeSignature}, ${forgedMessage}) is a valid signature!`);
console.log('');
console.log('  The attacker created a valid signature without the private key!');
console.log('  This is called "existential forgery".');
console.log('');
console.log('  Solution: Sign H(m) instead of m.');
console.log('  Attacker would need to find m such that H(m) = s^e mod n,');
console.log('  which requires inverting the hash function (infeasible).');

// ============================================================================
// Hash-and-Sign RSA
// ============================================================================

console.log('\n=== Hash-and-Sign RSA ===\n');

/**
 * Proper RSA signature:
 *   Sign:   s = H(m)^d mod n
 *   Verify: H(m) = s^e mod n
 *
 * The hash function H must be:
 * - Collision resistant: Hard to find m1 != m2 with H(m1) = H(m2)
 * - Preimage resistant: Hard to find m given H(m)
 * - Second preimage resistant: Given m1, hard to find m2 with H(m1) = H(m2)
 */

// Simple hash simulation (in practice, use SHA-256 or similar)
function simpleHash(data: string): bigint {
  // This is NOT a secure hash - just for demonstration
  let hash = 0n;
  for (let i = 0; i < data.length; i++) {
    hash = (hash * 31n + BigInt(data.charCodeAt(i))) % (2n ** 64n);
  }
  return hash;
}

const document = 'I owe Bob $100';
console.log(`Document: "${document}"`);

const hash = simpleHash(document);
console.log(`Hash: H("${document}") = ${hash}`);

// Reduce hash to fit within RSA modulus
const hashMod = hash % n;
console.log(`Hash mod n: ${hashMod}`);

// Sign the hash
const docSignature = Mod(hashMod, n).pow(d).value;
console.log(`\nSignature: s = H(m)^d mod n = ${docSignature}`);

// Verify
const recoveredHash = Mod(docSignature, n).pow(e).value;
console.log(`Verification: s^e mod n = ${recoveredHash}`);
console.log(`Matches original hash: ${recoveredHash === hashMod}`);

// ============================================================================
// Signature Forgery Prevention
// ============================================================================

console.log('\n=== Signature Forgery with Hashing ===\n');

console.log('Attempting existential forgery with hashing:');
console.log('');

const fakeS = 456n;
const forgedH = Mod(fakeS, n).pow(e).value;

console.log(`  1. Attacker chooses s = ${fakeS}`);
console.log(`  2. Attacker computes s^e mod n = ${forgedH}`);
console.log(`  3. For valid forgery, need document m where H(m) = ${forgedH}`);
console.log(`  4. Finding such m requires inverting hash function!`);
console.log(`  5. Attack fails: cannot find meaningful m with correct hash.`);

// ============================================================================
// Multiple Message Signing
// ============================================================================

console.log('\n=== Signing Multiple Messages ===\n');

const messages = ['Hello, World!', 'Transfer $1000 to Alice', 'Contract v1.0', ''];

console.log('Signing multiple documents:');
console.log('-'.repeat(50));

for (const msg of messages) {
  const h = simpleHash(msg) % n;
  const s = Mod(h, n).pow(d).value;
  const verified = Mod(s, n).pow(e).value === h;

  console.log(`Message: "${msg || '(empty)'}"`);
  console.log(`  Hash: ${h}`);
  console.log(`  Signature: ${s}`);
  console.log(`  Verified: ${verified}`);
  console.log('');
}

// ============================================================================
// Signature Malleability
// ============================================================================

console.log('\n=== Signature Malleability ===\n');

/**
 * RSA signatures can be malleable: given a valid signature,
 * an attacker can sometimes create another valid signature
 * for a related message.
 *
 * Example: Given signatures s1 for m1 and s2 for m2,
 * attacker can create signature s1*s2 for m1*m2.
 */

console.log('Multiplicative property of RSA signatures:');
console.log('');

const m1 = 10n;
const m2 = 20n;
const s1 = Mod(m1, n).pow(d).value;
const s2 = Mod(m2, n).pow(d).value;

console.log(`  m1 = ${m1}, s1 = m1^d = ${s1}`);
console.log(`  m2 = ${m2}, s2 = m2^d = ${s2}`);

// Compute s3 = s1 * s2 mod n
const s3 = (s1 * s2) % n;
const m3 = (m1 * m2) % n;

console.log(`\n  Combined signature: s3 = s1 * s2 mod n = ${s3}`);
console.log(`  Combined message: m3 = m1 * m2 mod n = ${m3}`);

// Verify s3 is valid signature for m3
const verified = Mod(s3, n).pow(e).value;
console.log(`\n  Verification: s3^e mod n = ${verified}`);
console.log(`  m3 = ${m3}`);
console.log(`  s3 is valid signature for m3: ${verified === m3}`);

console.log('\nAttack scenario:');
console.log('  Attacker has valid signatures for two "harmless" messages.');
console.log('  By multiplying them, attacker creates signature for m1*m2.');
console.log('  With careful choice, m1*m2 could be a "harmful" message!');

// ============================================================================
// RSA-PSS: Probabilistic Signature Scheme
// ============================================================================

console.log('\n=== RSA-PSS (Probabilistic Signature Scheme) ===\n');

/**
 * RSA-PSS is the recommended padding scheme for RSA signatures.
 * It adds randomness (salt) to make signatures non-deterministic.
 *
 * PSS Signature:
 * 1. Generate random salt s
 * 2. Compute m' = Hash(padding || Hash(M) || salt)
 * 3. Apply mask generation function
 * 4. Sign the padded value
 *
 * Benefits:
 * - Provably secure in the random oracle model
 * - Non-deterministic (different signature each time)
 * - Prevents various attacks on PKCS#1 v1.5
 */

console.log('RSA-PSS structure (conceptual):');
console.log('');
console.log('  Signing:');
console.log('    1. Generate random salt');
console.log('    2. DB = padding || salt');
console.log('    3. H = Hash(padding || Hash(message) || salt)');
console.log('    4. dbMask = MGF(H, len(DB))');
console.log('    5. maskedDB = DB XOR dbMask');
console.log('    6. EM = maskedDB || H || 0xbc');
console.log('    7. s = EM^d mod n');
console.log('');
console.log('  Verification:');
console.log('    1. Compute EM = s^e mod n');
console.log('    2. Parse EM to get maskedDB and H');
console.log('    3. dbMask = MGF(H, len(maskedDB))');
console.log('    4. DB = maskedDB XOR dbMask');
console.log('    5. Extract salt from DB');
console.log('    6. H\' = Hash(padding || Hash(message) || salt)');
console.log("    7. Verify H' = H");

// Demonstrate non-deterministic signing
console.log('\nNon-deterministic signing (simulated):');
const docToSign = 'Important document';
const docHash = simpleHash(docToSign);

for (let i = 0; i < 3; i++) {
  // Simulate adding random salt
  const salt = BigInt(Math.floor(Math.random() * 1000000));
  const saltedHash = (docHash * 1000000n + salt) % n;
  const sig = Mod(saltedHash, n).pow(d).value;

  console.log(`  Signing attempt ${i + 1}:`);
  console.log(`    Salt: ${salt}`);
  console.log(`    Signature: ${sig}`);
}
console.log('  Different signature each time, all valid!');

// ============================================================================
// Signature vs Encryption: Key Usage
// ============================================================================

console.log('\n=== Signature vs Encryption Key Usage ===\n');

/**
 * IMPORTANT: Never use the same key for both signing and encryption!
 *
 * Reason: Different security requirements and attack surfaces.
 *
 * Attack example:
 * 1. Attacker sends ciphertext c to victim
 * 2. Victim decrypts and "signs" the result: s = c^d mod n
 * 3. Attacker now has m = c^d mod n (the decryption!)
 *
 * Solution: Maintain separate key pairs for encryption and signing.
 */

console.log('Why separate keys for signing and encryption:');
console.log('');
console.log('  Attack if same key used:');
console.log('    1. Attacker creates "document" that is actually ciphertext c');
console.log('    2. Victim "signs" it: s = c^d mod n');
console.log('    3. But c^d is the decryption of c!');
console.log('    4. Attacker obtains decryption without access to private key');
console.log('');
console.log('  Best practice:');
console.log('    - One key pair for encryption: (n_enc, e_enc, d_enc)');
console.log('    - One key pair for signing: (n_sig, e_sig, d_sig)');
console.log('    - Certificate specifies key usage');

// ============================================================================
// RSA Blind Signatures
// ============================================================================

console.log('\n=== RSA Blind Signatures ===\n');

/**
 * Blind signatures allow someone to get a signature without
 * revealing what they're signing. Useful for:
 * - Anonymous digital cash
 * - Secret voting
 * - Credential issuance
 *
 * Protocol:
 * 1. User blinds message: m' = m * r^e mod n (r is random)
 * 2. Signer signs blinded message: s' = (m')^d mod n
 * 3. User unblinds: s = s' * r^(-1) mod n
 * 4. Result: s is valid signature on m!
 */

console.log('Blind signature protocol:');
console.log('');

const secretMessage = 777n;
console.log(`User's secret message: m = ${secretMessage}`);

// Step 1: User blinds the message
const blindingFactor = 123n; // Random r (must be coprime to n)
console.log(`Blinding factor: r = ${blindingFactor}`);

const rToE = Mod(blindingFactor, n).pow(e).value;
const blindedMessage = (secretMessage * rToE) % n;
console.log(`Blinded message: m' = m * r^e mod n = ${blindedMessage}`);

// Step 2: Signer signs the blinded message (without seeing real message)
const blindSignature = Mod(blindedMessage, n).pow(d).value;
console.log(`\nSigner produces: s' = (m')^d mod n = ${blindSignature}`);
console.log(`(Signer never sees the actual message ${secretMessage})`);

// Step 3: User unblinds the signature
const rInverse = Mod(blindingFactor, n).inv().value;
const unblindedSignature = (blindSignature * rInverse) % n;
console.log(`\nUser unblinds: s = s' * r^(-1) mod n = ${unblindedSignature}`);

// Verify: s should be valid signature on original message
const verification = Mod(unblindedSignature, n).pow(e).value;
console.log(`\nVerification: s^e mod n = ${verification}`);
console.log(`Original message: m = ${secretMessage}`);
console.log(`Signature valid: ${verification === secretMessage}`);

console.log('\nBlind signature properties:');
console.log('  - Signer cannot link signature to signing request');
console.log('  - Signature is valid for the unblinded message');
console.log('  - User proves knowledge of valid signature');

// ============================================================================
// Summary
// ============================================================================

console.log('\n=== Summary ===\n');

console.log('RSA SIGNATURES KEY POINTS:');
console.log('');
console.log('  1. Basic operation:');
console.log('     Sign:   s = H(m)^d mod n');
console.log('     Verify: H(m) = s^e mod n');
console.log('');
console.log('  2. Security requirements:');
console.log('     - Always sign hashes, not raw messages');
console.log('     - Use RSA-PSS padding scheme');
console.log('     - Use collision-resistant hash (SHA-256+)');
console.log('');
console.log('  3. Best practices:');
console.log('     - Separate keys for signing and encryption');
console.log('     - Minimum 2048-bit keys');
console.log('     - Consider ECDSA for better performance');
console.log('');
console.log('  4. Advanced applications:');
console.log('     - Blind signatures for anonymity');
console.log('     - Multi-signatures for threshold schemes');

console.log('\n' + '='.repeat(60));
console.log('RSA Digital Signatures Tutorial Complete');
console.log('='.repeat(60));
