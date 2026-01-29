/**
 * # ElGamal Encryption
 *
 * ElGamal encryption is a public-key encryption scheme based on the
 * Diffie-Hellman key exchange. It was proposed by Taher ElGamal in 1985.
 *
 * ## Key Insight
 *
 * ElGamal is "Diffie-Hellman with message blinding":
 * 1. Perform ephemeral DH to get shared secret K = g^(xr)
 * 2. Use K to blind the message: c2 = m * K
 *
 * ## Security
 *
 * - IND-CPA secure under DDH assumption
 * - Ciphertext is 2x message size (expansion factor 2)
 * - Randomized: encrypting same message gives different ciphertexts
 *
 * ## Properties
 *
 * - Homomorphic: E(m1) * E(m2) = E(m1 * m2)
 * - Probabilistic encryption
 * - Requires random number generation during encryption
 */

import { gcd, is_prime, power_mod, inverse_mod, factor, Mod } from 'sagemath-ts';

console.log('='.repeat(60));
console.log('ElGamal Encryption Tutorial');
console.log('='.repeat(60));

// ============================================================================
// Setup: ElGamal Key Generation
// ============================================================================

console.log('\n=== ElGamal Key Generation ===\n');

/**
 * ElGamal Key Generation:
 *
 * 1. Choose prime p and generator g
 * 2. Choose random private key x in [1, p-2]
 * 3. Compute public key y = g^x mod p
 *
 * Public key: (p, g, y)
 * Private key: x
 */

// For educational purposes, use small prime
const p = 23n;

// Find a generator
function findGenerator(p: bigint): bigint {
  const order = p - 1n;
  const factors = factor(order).filter(([f]) => f > 0n);

  for (let g = 2n; g < p; g++) {
    let isGenerator = true;

    for (const [q] of factors) {
      const exp = order / q;
      if (power_mod(g, exp, p) === 1n) {
        isGenerator = false;
        break;
      }
    }

    if (isGenerator) {
      return g;
    }
  }

  throw new Error('No generator found');
}

const g = findGenerator(p);

// Private key
const x = 6n; // Private key (random in [1, p-2])

// Public key
const y = Mod(g, p).pow(x).value;

console.log('Key Generation:');
console.log(`  Prime: p = ${p}`);
console.log(`  Generator: g = ${g}`);
console.log(`  Private key: x = ${x}`);
console.log(`  Public key: y = g^x = ${g}^${x} mod ${p} = ${y}`);

console.log('\nPublic Key (shared): (p, g, y) = (${p}, ${g}, ${y})');
console.log(`Private Key (secret): x = ${x}`);

// ============================================================================
// ElGamal Encryption
// ============================================================================

console.log('\n=== ElGamal Encryption ===\n');

/**
 * ElGamal Encryption:
 *
 * To encrypt message m (0 < m < p):
 * 1. Choose random r in [1, p-2]
 * 2. Compute c1 = g^r mod p
 * 3. Compute c2 = m * y^r mod p
 *
 * Ciphertext: (c1, c2)
 *
 * Understanding:
 * - c1 = g^r is ephemeral DH public value
 * - y^r = (g^x)^r = g^(xr) is shared DH secret
 * - c2 = m * g^(xr) blinds the message with shared secret
 */

const m = 7n; // Message to encrypt (must be in (Z/pZ)*)
console.log(`Message to encrypt: m = ${m}`);

// Encryption
const r = 3n; // Random ephemeral key (in practice, use secure random)
console.log(`Random ephemeral key: r = ${r}`);

const c1 = Mod(g, p).pow(r).value;
console.log(`\nStep 1: c1 = g^r mod p = ${g}^${r} mod ${p} = ${c1}`);

const y_r = Mod(y, p).pow(r).value;
console.log(`Step 2: Compute shared secret y^r = ${y}^${r} mod ${p} = ${y_r}`);

const c2 = (m * y_r) % p;
console.log(`Step 3: c2 = m * y^r mod p = ${m} * ${y_r} mod ${p} = ${c2}`);

console.log(`\nCiphertext: (c1, c2) = (${c1}, ${c2})`);

// ============================================================================
// ElGamal Decryption
// ============================================================================

console.log('\n=== ElGamal Decryption ===\n');

/**
 * ElGamal Decryption:
 *
 * Given ciphertext (c1, c2) and private key x:
 * 1. Compute shared secret s = c1^x mod p
 * 2. Compute message m = c2 * s^(-1) mod p
 *
 * Why this works:
 * - c1 = g^r, so c1^x = g^(rx) = g^(xr)
 * - c2 = m * g^(xr), so c2 * (g^(xr))^(-1) = m
 */

console.log('Decryption:');
console.log(`  Ciphertext: (c1, c2) = (${c1}, ${c2})`);
console.log(`  Private key: x = ${x}`);

// Step 1: Compute shared secret
const s = Mod(c1, p).pow(x).value;
console.log(`\nStep 1: s = c1^x mod p = ${c1}^${x} mod ${p} = ${s}`);

// Step 2: Compute inverse of shared secret
const s_inv = Mod(s, p).inv().value;
console.log(`Step 2: s^(-1) mod p = ${s_inv}`);

// Step 3: Recover message
const m_recovered = (c2 * s_inv) % p;
console.log(`Step 3: m = c2 * s^(-1) mod p = ${c2} * ${s_inv} mod ${p} = ${m_recovered}`);

console.log(`\nOriginal message: ${m}`);
console.log(`Recovered message: ${m_recovered}`);
console.log(`Decryption successful: ${m === m_recovered}`);

// ============================================================================
// Why ElGamal Works
// ============================================================================

console.log('\n=== Why ElGamal Works ===\n');

console.log('Mathematical justification:');
console.log('');
console.log('  Encryption:');
console.log(`    c1 = g^r = ${c1}`);
console.log(`    c2 = m * y^r = m * (g^x)^r = m * g^(xr) = ${c2}`);
console.log('');
console.log('  Decryption:');
console.log(`    s = c1^x = (g^r)^x = g^(rx) = g^(xr) = ${s}`);
console.log(`    m = c2 * s^(-1) = (m * g^(xr)) * (g^(xr))^(-1) = m = ${m_recovered}`);
console.log('');
console.log('  Key observation: Both sides compute g^(xr) as shared secret.');
console.log('  This is exactly the DH shared secret!');

// ============================================================================
// Probabilistic Encryption
// ============================================================================

console.log('\n=== Probabilistic Encryption ===\n');

/**
 * ElGamal is probabilistic: same message encrypts to different ciphertexts.
 * This is essential for semantic security (IND-CPA).
 */

console.log(`Encrypting m = ${m} with different random values:`);
console.log('');

for (let r_test = 1n; r_test <= 5n; r_test++) {
  const c1_test = Mod(g, p).pow(r_test).value;
  const c2_test = (m * Mod(y, p).pow(r_test).value) % p;
  console.log(`  r = ${r_test}: (c1, c2) = (${c1_test}, ${c2_test})`);
}

console.log('');
console.log('Each encryption is different, even for the same message!');
console.log('This prevents chosen-plaintext attacks.');

// ============================================================================
// Homomorphic Property
// ============================================================================

console.log('\n=== Homomorphic Property ===\n');

/**
 * ElGamal is multiplicatively homomorphic:
 *   E(m1) * E(m2) = E(m1 * m2)
 *
 * Given (c1, c2) = E(m1) and (c1', c2') = E(m2):
 *   (c1 * c1', c2 * c2') = E(m1 * m2)
 *
 * This enables computation on encrypted data!
 */

console.log('Homomorphic multiplication:');
console.log('');

// Encrypt two messages
const m1 = 3n;
const m2 = 4n;
const r1 = 5n;
const r2 = 7n;

const c1_m1 = Mod(g, p).pow(r1).value;
const c2_m1 = (m1 * Mod(y, p).pow(r1).value) % p;

const c1_m2 = Mod(g, p).pow(r2).value;
const c2_m2 = (m2 * Mod(y, p).pow(r2).value) % p;

console.log(`E(m1=${m1}) = (${c1_m1}, ${c2_m1})`);
console.log(`E(m2=${m2}) = (${c1_m2}, ${c2_m2})`);

// Multiply ciphertexts
const c1_product = (c1_m1 * c1_m2) % p;
const c2_product = (c2_m1 * c2_m2) % p;

console.log(`\nE(m1) * E(m2) = (${c1_product}, ${c2_product})`);

// Decrypt the product
const s_product = Mod(c1_product, p).pow(x).value;
const m_product = (c2_product * Mod(s_product, p).inv().value) % p;

console.log(`Decryption: ${m_product}`);
console.log(`Expected m1 * m2 = ${m1} * ${m2} = ${(m1 * m2) % p}`);
console.log(`Homomorphism verified: ${m_product === (m1 * m2) % p}`);

console.log('\nImplications:');
console.log('  - Can compute on encrypted data');
console.log('  - Used in e-voting, secure auctions');
console.log('  - Foundation for more complex homomorphic schemes');

// ============================================================================
// Security Analysis
// ============================================================================

console.log('\n=== Security Analysis ===\n');

/**
 * ElGamal security under DDH:
 *
 * If DDH is hard, then ElGamal is IND-CPA secure.
 *
 * Proof sketch:
 * - Real: (g, g^x, g^r, g^(xr))
 * - Fake: (g, g^x, g^r, g^z) for random z
 * - Ciphertext uses g^(xr), which is indistinguishable from random
 * - Therefore, c2 = m * g^(xr) hides m
 */

console.log('Security under DDH assumption:');
console.log('');
console.log('  If adversary can distinguish:');
console.log('    (g, y=g^x, c1=g^r, c2=m*g^(xr))');
console.log('  from');
console.log('    (g, y=g^x, c1=g^r, c2=m*g^z)  [random z]');
console.log('');
console.log('  Then adversary can distinguish DDH tuples:');
console.log('    (g, g^x, g^r, g^(xr)) from (g, g^x, g^r, g^z)');
console.log('');
console.log('  Since DDH is hard, ElGamal is IND-CPA secure.');

// ============================================================================
// Malleability Warning
// ============================================================================

console.log('\n=== Malleability Warning ===\n');

/**
 * ElGamal is malleable: given E(m), attacker can create E(k*m) for known k.
 *
 * Attack: Given (c1, c2) = E(m), compute (c1, c2 * k) = E(m * k)
 *
 * This means ElGamal is NOT CCA-secure (chosen ciphertext attack).
 * Use hybrid encryption or ElGamal-KEM for CCA security.
 */

console.log('Malleability attack:');
console.log('');

const k_attack = 5n;
console.log(`  Original ciphertext E(${m}) = (${c1}, ${c2})`);

// Attacker modifies ciphertext
const c2_modified = (c2 * k_attack) % p;
console.log(`  Attacker multiplies c2 by k=${k_attack}`);
console.log(`  Modified ciphertext: (${c1}, ${c2_modified})`);

// Decrypt modified ciphertext
const s_modified = Mod(c1, p).pow(x).value;
const m_modified = (c2_modified * Mod(s_modified, p).inv().value) % p;

console.log(`\n  Decryption of modified ciphertext: ${m_modified}`);
console.log(`  Expected: m * k = ${m} * ${k_attack} = ${(m * k_attack) % p}`);
console.log(`  Malleability demonstrated: ${m_modified === (m * k_attack) % p}`);

console.log('\nImplications:');
console.log('  - Attacker can change encrypted vote from 1 to 5');
console.log('  - Attacker can manipulate encrypted financial data');
console.log('  - Solution: Use authenticated encryption (ElGamal + MAC)');

// ============================================================================
// Complete ElGamal Class
// ============================================================================

console.log('\n=== ElGamal Implementation ===\n');

interface ElGamalPublicKey {
  p: bigint;
  g: bigint;
  y: bigint;
}

interface ElGamalPrivateKey {
  p: bigint;
  x: bigint;
}

interface ElGamalCiphertext {
  c1: bigint;
  c2: bigint;
}

function elgamalEncrypt(
  m: bigint,
  publicKey: ElGamalPublicKey,
  r?: bigint
): ElGamalCiphertext {
  const { p, g, y } = publicKey;

  // Generate random r if not provided
  if (r === undefined) {
    r = BigInt(Math.floor(Math.random() * (Number(p) - 2))) + 1n;
  }

  const c1 = power_mod(g, r, p);
  const c2 = (m * power_mod(y, r, p)) % p;

  return { c1, c2 };
}

function elgamalDecrypt(
  ciphertext: ElGamalCiphertext,
  privateKey: ElGamalPrivateKey
): bigint {
  const { c1, c2 } = ciphertext;
  const { p, x } = privateKey;

  const s = power_mod(c1, x, p);
  const s_inv = inverse_mod(s, p);

  return (c2 * s_inv) % p;
}

// Test the implementation
const pubKey: ElGamalPublicKey = { p, g, y };
const privKey: ElGamalPrivateKey = { p, x };

console.log('Testing ElGamal implementation:');

for (const testMsg of [1n, 5n, 10n, 15n, 20n]) {
  if (testMsg >= p) continue;

  const encrypted = elgamalEncrypt(testMsg, pubKey);
  const decrypted = elgamalDecrypt(encrypted, privKey);

  console.log(
    `  m=${testMsg}: E(m)=(${encrypted.c1}, ${encrypted.c2}), D(E(m))=${decrypted}, correct=${testMsg === decrypted}`
  );
}

// ============================================================================
// Encoding Messages
// ============================================================================

console.log('\n=== Encoding Messages ===\n');

/**
 * ElGamal encrypts elements of (Z/pZ)*, not arbitrary byte strings.
 * We need encoding schemes:
 *
 * 1. Simple: Interpret bytes as integer < p
 * 2. With padding: Add structure for verification
 * 3. Hybrid: Encrypt symmetric key, use symmetric cipher for data
 */

console.log('Message encoding considerations:');
console.log('');
console.log('  1. Message must be in (Z/pZ)*');
console.log('     - 0 < m < p');
console.log('     - gcd(m, p) = 1 (automatic for prime p if m > 0)');
console.log('');
console.log('  2. Encoding methods:');
console.log('     - Bytes to integer (for small messages)');
console.log('     - Hybrid encryption (for large messages)');
console.log('');
console.log('  3. Ciphertext expansion:');
console.log('     - Message: |m| bits');
console.log('     - Ciphertext: 2|p| bits (c1 and c2)');
console.log('     - Expansion factor: ~2x');

// ============================================================================
// Summary
// ============================================================================

console.log('\n=== Summary ===\n');

console.log('ELGAMAL ENCRYPTION KEY POINTS:');
console.log('');
console.log('  1. Key Generation:');
console.log('     - Choose prime p, generator g');
console.log('     - Private key: random x');
console.log('     - Public key: y = g^x');
console.log('');
console.log('  2. Encryption E(m):');
console.log('     - Choose random r');
console.log('     - c1 = g^r, c2 = m * y^r');
console.log('     - Ciphertext: (c1, c2)');
console.log('');
console.log('  3. Decryption:');
console.log('     - s = c1^x');
console.log('     - m = c2 * s^(-1)');
console.log('');
console.log('  4. Properties:');
console.log('     - Probabilistic (IND-CPA secure under DDH)');
console.log('     - Multiplicatively homomorphic');
console.log('     - Malleable (NOT CCA secure without modification)');
console.log('');
console.log('  5. Best Practices:');
console.log('     - Use hybrid encryption for long messages');
console.log('     - Add authentication for CCA security');
console.log('     - Consider ECIES for better efficiency');

console.log('\n' + '='.repeat(60));
console.log('ElGamal Encryption Tutorial Complete');
console.log('='.repeat(60));
