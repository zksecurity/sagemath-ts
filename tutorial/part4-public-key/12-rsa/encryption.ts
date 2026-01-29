/**
 * # RSA Encryption
 *
 * This tutorial covers RSA encryption: how to encrypt messages using the public key.
 *
 * ## The Encryption Function
 *
 * Given a public key (n, e) and a message m where 0 <= m < n:
 *   c = m^e mod n
 *
 * The ciphertext c can only be decrypted by someone who knows the private key d.
 *
 * ## Important Constraints
 *
 * 1. The message m must be less than the modulus n
 * 2. The message should be coprime to n (almost always true for random m)
 * 3. Raw RSA is deterministic - same message gives same ciphertext
 * 4. Padding schemes (OAEP, PKCS#1) are essential for security
 *
 * ## Security Considerations
 *
 * - Never use "textbook RSA" in practice
 * - Always use proper padding (OAEP for encryption)
 * - RSA is slow compared to symmetric encryption
 * - Hybrid encryption: Use RSA to encrypt a symmetric key
 */

import { gcd, inverse_mod, power_mod, Mod } from 'sagemath-ts';

console.log('='.repeat(60));
console.log('RSA Encryption Tutorial');
console.log('='.repeat(60));

// ============================================================================
// Setup: RSA Key Pair
// ============================================================================

console.log('\n=== RSA Key Setup ===\n');

// Key parameters (from key-generation.ts)
const p = 61n;
const q = 53n;
const n = p * q; // 3233
const phi_n = (p - 1n) * (q - 1n); // 3120
const e = 17n;
const d = inverse_mod(e, phi_n); // 2753

console.log('Public Key (shared with everyone):');
console.log(`  n = ${n}`);
console.log(`  e = ${e}`);

console.log('\nPrivate Key (kept secret):');
console.log(`  d = ${d}`);

// ============================================================================
// Basic RSA Encryption
// ============================================================================

console.log('\n=== Basic RSA Encryption ===\n');

/**
 * RSA Encryption:
 *   c = m^e mod n
 *
 * Where:
 * - m is the plaintext message (as an integer)
 * - e is the public exponent
 * - n is the modulus
 * - c is the ciphertext
 */

const message = 65n; // 'A' in ASCII
console.log(`Original message: m = ${message}`);
console.log(`(This represents ASCII 'A')`);

// Encrypt using modular exponentiation
const ciphertext = Mod(message, n).pow(e).value;
console.log(`\nEncryption: c = m^e mod n`);
console.log(`  c = ${message}^${e} mod ${n}`);
console.log(`  c = ${ciphertext}`);

// Verify decryption works
const decrypted = Mod(ciphertext, n).pow(d).value;
console.log(`\nDecryption verification: m' = c^d mod n`);
console.log(`  m' = ${ciphertext}^${d} mod ${n}`);
console.log(`  m' = ${decrypted}`);
console.log(`  Original recovered: ${decrypted === message}`);

// ============================================================================
// Encrypting Text Messages
// ============================================================================

console.log('\n=== Encrypting Text Messages ===\n');

/**
 * To encrypt text:
 * 1. Convert text to numbers (e.g., ASCII)
 * 2. Ensure each number < n (may need to break into blocks)
 * 3. Encrypt each number
 *
 * NOTE: This is for educational purposes only.
 * Real systems use proper padding schemes.
 */

function textToNumber(text: string): bigint {
  // Simple conversion: treat as base-256 number
  let result = 0n;
  for (let i = 0; i < text.length; i++) {
    result = result * 256n + BigInt(text.charCodeAt(i));
  }
  return result;
}

function numberToText(num: bigint): string {
  const chars: string[] = [];
  while (num > 0n) {
    chars.unshift(String.fromCharCode(Number(num % 256n)));
    num = num / 256n;
  }
  return chars.join('');
}

// Encrypt a short message (must be < n)
const shortMessage = 'HI';
const messageNum = textToNumber(shortMessage);

console.log(`Text message: "${shortMessage}"`);
console.log(`As number: ${messageNum}`);

if (messageNum >= n) {
  console.log(`ERROR: Message too large for this modulus!`);
} else {
  const encryptedNum = Mod(messageNum, n).pow(e).value;
  console.log(`Encrypted: ${encryptedNum}`);

  const decryptedNum = Mod(encryptedNum, n).pow(d).value;
  const decryptedText = numberToText(decryptedNum);
  console.log(`Decrypted number: ${decryptedNum}`);
  console.log(`Decrypted text: "${decryptedText}"`);
}

// ============================================================================
// Block Encryption for Longer Messages
// ============================================================================

console.log('\n=== Block Encryption for Longer Messages ===\n');

/**
 * When the message is larger than n, we must split it into blocks.
 * Each block is encrypted independently.
 *
 * Block size = floor(log2(n)) bits
 * For n = 3233, block size = 11 bits (but we'll use byte-aligned blocks)
 */

// Calculate max bytes per block
const maxValue = n - 1n;
let maxBytes = 0;
let temp = maxValue;
while (temp > 0n) {
  temp = temp / 256n;
  maxBytes++;
}
maxBytes--; // Ensure the block value is strictly < n

console.log(`Modulus n = ${n}`);
console.log(`Max bytes per block: ${maxBytes}`);

// For our small n, we can only fit 1 byte per block
function encryptBlocks(plaintext: string, n: bigint, e: bigint, blockSize: number): bigint[] {
  const cipherBlocks: bigint[] = [];

  for (let i = 0; i < plaintext.length; i += blockSize) {
    const block = plaintext.slice(i, i + blockSize);
    const blockNum = textToNumber(block);

    if (blockNum >= n) {
      throw new Error(`Block "${block}" (${blockNum}) >= n (${n})`);
    }

    const encryptedBlock = Mod(blockNum, n).pow(e).value;
    cipherBlocks.push(encryptedBlock);
  }

  return cipherBlocks;
}

function decryptBlocks(cipherBlocks: bigint[], n: bigint, d: bigint): string {
  let plaintext = '';

  for (const block of cipherBlocks) {
    const decryptedNum = Mod(block, n).pow(d).value;
    plaintext += numberToText(decryptedNum);
  }

  return plaintext;
}

const longMessage = 'HELLO';
console.log(`\nMessage: "${longMessage}"`);

const blocks = encryptBlocks(longMessage, n, e, 1);
console.log(`Encrypted blocks: [${blocks.join(', ')}]`);

const recovered = decryptBlocks(blocks, n, d);
console.log(`Decrypted: "${recovered}"`);

// ============================================================================
// Why Textbook RSA is Insecure
// ============================================================================

console.log('\n=== Why Textbook RSA is Insecure ===\n');

/**
 * "Textbook RSA" (RSA without padding) has several weaknesses:
 *
 * 1. DETERMINISTIC: Same plaintext always gives same ciphertext
 *    - Attacker can build a dictionary of common messages
 *
 * 2. MALLEABILITY: Given c = m^e mod n, attacker can compute:
 *    - c' = c * 2^e mod n = (2m)^e mod n
 *    - Without knowing m, they created valid encryption of 2m
 *
 * 3. SMALL MESSAGE ATTACK: If m^e < n, then c = m^e (no modular reduction)
 *    - Attacker can just take e-th root to recover m
 */

console.log('1. DETERMINISM PROBLEM:');
const msg1 = 100n;
const cipher1a = Mod(msg1, n).pow(e).value;
const cipher1b = Mod(msg1, n).pow(e).value;
console.log(`   Encrypting ${msg1} twice:`);
console.log(`   First encryption:  ${cipher1a}`);
console.log(`   Second encryption: ${cipher1b}`);
console.log(`   Same ciphertext: ${cipher1a === cipher1b}`);

console.log('\n2. MALLEABILITY PROBLEM:');
const originalMsg = 50n;
const originalCipher = Mod(originalMsg, n).pow(e).value;
console.log(`   Original: m = ${originalMsg}, c = ${originalCipher}`);

// Attacker modifies ciphertext without knowing m
const factor = 2n;
const factorEncrypted = Mod(factor, n).pow(e).value;
const modifiedCipher = (originalCipher * factorEncrypted) % n;
console.log(`   Attacker computes: c' = c * ${factor}^e mod n = ${modifiedCipher}`);

// When decrypted, this gives 2*m
const modifiedDecrypted = Mod(modifiedCipher, n).pow(d).value;
console.log(`   Decrypted c': ${modifiedDecrypted}`);
console.log(`   This equals ${factor} * ${originalMsg} = ${factor * originalMsg}!`);

console.log('\n3. SMALL MESSAGE WITH SMALL EXPONENT:');
// With e=3 and a small message, m^e might not wrap around
const smallE = 3n;
const smallMsg = 10n;
const rawPower = smallMsg ** smallE;
const smallN = 3233n;
console.log(`   If e = ${smallE} and m = ${smallMsg}:`);
console.log(`   m^e = ${rawPower}`);
console.log(`   With n = ${smallN}, c = ${rawPower % smallN}`);
if (rawPower < smallN) {
  console.log(`   Since m^e = ${rawPower} < n = ${smallN}, the ciphertext IS m^e`);
  console.log(`   Attacker can compute cube root to get m!`);
}

// ============================================================================
// OAEP Padding (Conceptual Overview)
// ============================================================================

console.log('\n=== OAEP Padding (Conceptual) ===\n');

/**
 * OAEP (Optimal Asymmetric Encryption Padding) adds randomness and structure:
 *
 * 1. Pad message m with zeros to fixed length
 * 2. Generate random r
 * 3. Compute: X = (m || zeros) XOR Hash(r)
 * 4. Compute: Y = r XOR Hash(X)
 * 5. Encrypt (X || Y) with RSA
 *
 * Benefits:
 * - Randomness makes encryption non-deterministic
 * - Structure allows detection of tampering
 * - Proven secure under random oracle model
 */

console.log('OAEP structure:');
console.log('  1. Generate random seed r');
console.log('  2. X = (m || padding) XOR G(r)     // G is hash function');
console.log('  3. Y = r XOR H(X)                  // H is another hash function');
console.log('  4. Encrypt (X || Y)');

console.log('\nDecryption reverses the process:');
console.log('  1. Decrypt to get (X || Y)');
console.log('  2. r = Y XOR H(X)');
console.log('  3. m || padding = X XOR G(r)');
console.log('  4. Verify and remove padding');

// Demonstrate that adding randomness makes encryption non-deterministic
console.log('\nSimulated non-deterministic encryption:');
for (let i = 0; i < 3; i++) {
  // Simulate adding random padding
  const randomPadding = BigInt(Math.floor(Math.random() * 1000));
  const paddedMsg = (50n * 1000n + randomPadding) % n;
  const cipher = Mod(paddedMsg, n).pow(e).value;
  console.log(`  Encryption ${i + 1}: padding=${randomPadding}, ciphertext=${cipher}`);
}

// ============================================================================
// Hybrid Encryption
// ============================================================================

console.log('\n=== Hybrid Encryption ===\n');

/**
 * RSA is slow (complex modular exponentiation) and has message size limits.
 * In practice, we use HYBRID ENCRYPTION:
 *
 * 1. Generate a random symmetric key K
 * 2. Encrypt message with symmetric cipher: C_sym = AES(K, message)
 * 3. Encrypt key with RSA: C_key = RSA(K)
 * 4. Send (C_key, C_sym)
 *
 * Benefits:
 * - RSA only encrypts small key (fast)
 * - AES encrypts large message (fast)
 * - Combines benefits of both systems
 */

console.log('Hybrid Encryption Process:');
console.log('');
console.log('  Sender:');
console.log('    1. Generate random 256-bit symmetric key K');
console.log('    2. C_message = AES-GCM(K, plaintext)');
console.log('    3. C_key = RSA-OAEP(public_key, K)');
console.log('    4. Send (C_key, C_message)');
console.log('');
console.log('  Receiver:');
console.log('    1. K = RSA-OAEP-decrypt(private_key, C_key)');
console.log('    2. plaintext = AES-GCM-decrypt(K, C_message)');

// Simulate hybrid encryption
console.log('\nSimulation with our RSA:');
const symmetricKey = 2048n; // Simulated AES key (would be 256 bits)
console.log(`  Symmetric key K = ${symmetricKey}`);

const encryptedKey = Mod(symmetricKey, n).pow(e).value;
console.log(`  RSA-encrypted key = ${encryptedKey}`);

const decryptedKey = Mod(encryptedKey, n).pow(d).value;
console.log(`  Decrypted key = ${decryptedKey}`);
console.log(`  Key recovered successfully: ${decryptedKey === symmetricKey}`);

// ============================================================================
// Encryption Performance Considerations
// ============================================================================

console.log('\n=== Performance Considerations ===\n');

/**
 * RSA encryption speed depends on:
 * 1. Size of the exponent e
 * 2. Size of the modulus n
 * 3. Implementation (square-and-multiply optimization)
 *
 * Using e = 65537 (binary: 10000000000000001):
 * - Only 2 set bits (Hamming weight 2)
 * - Square-and-multiply requires only ~17 squarings and 2 multiplications
 * - Much faster than using a large random e
 */

console.log('Why e = 65537 is fast:');
const e65537 = 65537n;
console.log(`  e = ${e65537}`);
console.log(`  Binary: ${e65537.toString(2)}`);
console.log(`  Bit length: ${e65537.toString(2).length}`);
console.log(`  Number of 1-bits: ${e65537.toString(2).split('').filter((b) => b === '1').length}`);
console.log('');
console.log('  Square-and-multiply for m^65537:');
console.log('    - 16 squarings to compute m^(2^16)');
console.log('    - 1 multiplication: m^(2^16) * m = m^65537');
console.log('    Total: 17 modular operations');

// ============================================================================
// Summary
// ============================================================================

console.log('\n=== Summary ===\n');

console.log('RSA ENCRYPTION KEY POINTS:');
console.log('');
console.log('  1. Basic operation: c = m^e mod n');
console.log('');
console.log('  2. Constraints:');
console.log('     - Message must be < n');
console.log('     - Larger messages need block encryption');
console.log('');
console.log('  3. Security requirements:');
console.log('     - NEVER use textbook RSA');
console.log('     - Always use OAEP or similar padding');
console.log('     - Use hybrid encryption for efficiency');
console.log('');
console.log('  4. Standard practice:');
console.log('     - e = 65537 for fast encryption');
console.log('     - RSA-OAEP for padding');
console.log('     - Encrypt symmetric keys, not data');

console.log('\n' + '='.repeat(60));
console.log('RSA Encryption Tutorial Complete');
console.log('='.repeat(60));
