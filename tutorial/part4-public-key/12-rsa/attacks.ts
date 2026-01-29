/**
 * # RSA Attacks
 *
 * This tutorial demonstrates common attacks on RSA implementations and
 * explains why certain parameter choices or practices are dangerous.
 *
 * ## Attack Categories
 *
 * 1. Mathematical attacks (small e, common modulus, etc.)
 * 2. Implementation attacks (timing, padding oracle)
 * 3. Key generation attacks (weak primes, close primes)
 *
 * Understanding these attacks is crucial for secure RSA implementation.
 *
 * ## DISCLAIMER
 *
 * These attacks are presented for educational purposes only.
 * Use this knowledge to build more secure systems, not to attack others.
 */

import { gcd, inverse_mod, isqrt, factor, xgcd, crt, Mod } from 'sagemath-ts';

console.log('='.repeat(60));
console.log('RSA Attacks Tutorial');
console.log('='.repeat(60));

// ============================================================================
// Attack 1: Small Public Exponent with Small Message
// ============================================================================

console.log('\n=== Attack 1: Small Public Exponent Attack ===\n');

/**
 * If e is small (e.g., e=3) and the message m is small enough that
 * m^e < n, then the ciphertext c = m^e without any modular reduction.
 *
 * Attack: Simply compute the e-th root of c.
 */

console.log('Scenario: e = 3, small message');
console.log('');

const n1 = 3233n; // Our example modulus
const e_small = 3n;

// A message small enough that m^3 < n
const smallMessage = 10n;
const smallCubed = smallMessage ** e_small;

console.log(`n = ${n1}`);
console.log(`e = ${e_small}`);
console.log(`m = ${smallMessage}`);
console.log(`m^e = ${smallCubed}`);
console.log(`m^e < n? ${smallCubed < n1}`);

if (smallCubed < n1) {
  // Ciphertext is just m^e, no modular reduction happened
  const ciphertext = smallCubed;
  console.log(`\nCiphertext c = ${ciphertext} (no modular reduction)`);

  // Attack: compute cube root
  // For educational purposes, we compute integer cube root
  function integerNthRoot(x: bigint, n: number): bigint {
    if (x < 0n) throw new Error('x must be non-negative');
    if (x === 0n) return 0n;
    if (n === 1) return x;

    let low = 1n;
    let high = x;

    while (low < high) {
      const mid = (low + high + 1n) / 2n;
      if (mid ** BigInt(n) <= x) {
        low = mid;
      } else {
        high = mid - 1n;
      }
    }

    return low;
  }

  const recovered = integerNthRoot(ciphertext, 3);
  console.log(`Attack: cube root of ${ciphertext} = ${recovered}`);
  console.log(`Original message recovered: ${recovered === smallMessage}`);
}

console.log('\nMitigation:');
console.log('  - Use proper padding (OAEP) which adds randomness');
console.log('  - Use e = 65537 instead of e = 3');
console.log('  - Ensure padded message is always close to n in size');

// ============================================================================
// Attack 2: Hastad's Broadcast Attack
// ============================================================================

console.log('\n=== Attack 2: Hastad\'s Broadcast Attack ===\n');

/**
 * If the same message m is encrypted with e different public keys,
 * all using exponent e, the message can be recovered using CRT.
 *
 * Given:
 *   c1 = m^e mod n1
 *   c2 = m^e mod n2
 *   c3 = m^e mod n3   (for e=3)
 *
 * Find x such that:
 *   x ≡ c1 (mod n1)
 *   x ≡ c2 (mod n2)
 *   x ≡ c3 (mod n3)
 *
 * By CRT, x = m^e (no modular reduction if m^e < n1*n2*n3).
 * Then m = e-th root of x.
 */

console.log('Scenario: Same message sent to 3 recipients with e=3');
console.log('');

// Three different RSA moduli (pairwise coprime)
const n_1 = 3233n; // 61 * 53
const n_2 = 2773n; // 59 * 47
const n_3 = 3127n; // 53 * 59

const e_broadcast = 3n;
const secretMessage = 10n; // Must be small enough that m^3 < n1*n2*n3

console.log('Recipients:');
console.log(`  n1 = ${n_1}`);
console.log(`  n2 = ${n_2}`);
console.log(`  n3 = ${n_3}`);
console.log(`  All use e = ${e_broadcast}`);
console.log(`\nSecret message: m = ${secretMessage}`);

// Encrypt for each recipient
const c_1 = Mod(secretMessage, n_1).pow(e_broadcast).value;
const c_2 = Mod(secretMessage, n_2).pow(e_broadcast).value;
const c_3 = Mod(secretMessage, n_3).pow(e_broadcast).value;

console.log('\nCiphertexts intercepted by attacker:');
console.log(`  c1 = ${c_1}`);
console.log(`  c2 = ${c_2}`);
console.log(`  c3 = ${c_3}`);

// Attacker uses CRT to find m^3
function CRT_list_simple(residues: bigint[], moduli: bigint[]): bigint {
  let result = residues[0]!;
  let modulus = moduli[0]!;

  for (let i = 1; i < residues.length; i++) {
    result = crt(result, residues[i]!, modulus, moduli[i]!);
    const [g] = xgcd(modulus, moduli[i]!);
    modulus = (modulus * moduli[i]!) / g;
  }

  return result;
}

const combined = CRT_list_simple([c_1, c_2, c_3], [n_1, n_2, n_3]);
console.log(`\nCRT combination: m^3 = ${combined}`);

// Check if we can take the cube root
const N_product = n_1 * n_2 * n_3;
const m_cubed = secretMessage ** 3n;
console.log(`m^3 = ${m_cubed}`);
console.log(`n1 * n2 * n3 = ${N_product}`);
console.log(`m^3 < n1*n2*n3? ${m_cubed < N_product}`);

// Take cube root to recover message
function integerCubeRoot(x: bigint): bigint {
  if (x < 0n) throw new Error('x must be non-negative');
  if (x < 2n) return x;

  let low = 1n;
  let high = x;

  while (low < high) {
    const mid = (low + high + 1n) / 2n;
    if (mid * mid * mid <= x) {
      low = mid;
    } else {
      high = mid - 1n;
    }
  }

  return low;
}

const attackerRecovered = integerCubeRoot(combined);
console.log(`\nAttacker computes cube root: ${attackerRecovered}`);
console.log(`Message recovered: ${attackerRecovered === secretMessage}`);

console.log('\nMitigation:');
console.log('  - Use unique random padding for each encryption');
console.log('  - Use larger e (e.g., 65537)');
console.log('  - OAEP padding prevents this attack');

// ============================================================================
// Attack 3: Common Modulus Attack
// ============================================================================

console.log('\n=== Attack 3: Common Modulus Attack ===\n');

/**
 * If the same message is encrypted under two different public exponents
 * e1 and e2 with gcd(e1, e2) = 1, using the same modulus n,
 * the message can be recovered without the private key.
 *
 * Given c1 = m^e1 mod n and c2 = m^e2 mod n:
 * 1. Find s, t such that s*e1 + t*e2 = 1 (by extended GCD)
 * 2. Compute: c1^s * c2^t = m^(s*e1 + t*e2) = m^1 = m (mod n)
 */

console.log('Scenario: Same n, different e values');
console.log('');

const n_common = 3233n;
const e1 = 17n;
const e2 = 19n; // Different exponent, gcd(17, 19) = 1

console.log(`Common modulus: n = ${n_common}`);
console.log(`Exponent 1: e1 = ${e1}`);
console.log(`Exponent 2: e2 = ${e2}`);
console.log(`gcd(e1, e2) = ${gcd(e1, e2)}`);

const secret = 42n;
console.log(`\nSecret message: m = ${secret}`);

const c1_common = Mod(secret, n_common).pow(e1).value;
const c2_common = Mod(secret, n_common).pow(e2).value;

console.log(`Ciphertext 1: c1 = m^e1 mod n = ${c1_common}`);
console.log(`Ciphertext 2: c2 = m^e2 mod n = ${c2_common}`);

// Find s, t such that s*e1 + t*e2 = 1
const [g, s, t] = xgcd(e1, e2);
console.log(`\nExtended GCD: ${s}*${e1} + ${t}*${e2} = ${g}`);

// Compute m = c1^s * c2^t mod n
// Note: if s or t is negative, we use modular inverse
let result_common: bigint;

if (s >= 0n && t >= 0n) {
  result_common = (Mod(c1_common, n_common).pow(s).value * Mod(c2_common, n_common).pow(t).value) % n_common;
} else if (s < 0n) {
  // c1^s = (c1^(-1))^(-s)
  const c1_inv = Mod(c1_common, n_common).inv().value;
  result_common = (Mod(c1_inv, n_common).pow(-s).value * Mod(c2_common, n_common).pow(t).value) % n_common;
} else {
  // t < 0n
  const c2_inv = Mod(c2_common, n_common).inv().value;
  result_common = (Mod(c1_common, n_common).pow(s).value * Mod(c2_inv, n_common).pow(-t).value) % n_common;
}

console.log(`\nAttack: m = c1^s * c2^t mod n`);
console.log(`Recovered message: ${result_common}`);
console.log(`Correct: ${result_common === secret}`);

console.log('\nMitigation:');
console.log('  - Never share the same n between different key pairs');
console.log('  - Each RSA key should have unique n');

// ============================================================================
// Attack 4: Fermat Factorization (Close Primes)
// ============================================================================

console.log('\n=== Attack 4: Fermat Factorization ===\n');

/**
 * If p and q are close together, n can be factored quickly using
 * Fermat's factorization method.
 *
 * For n = p * q where p > q:
 *   n = ((p+q)/2)^2 - ((p-q)/2)^2
 *   n = a^2 - b^2 = (a+b)(a-b)
 *
 * Algorithm:
 *   a = ceil(sqrt(n))
 *   while a^2 - n is not a perfect square:
 *     a += 1
 *   b = sqrt(a^2 - n)
 *   p = a + b, q = a - b
 *
 * If p and q are close, this converges very quickly.
 */

console.log('Scenario: Primes p and q are close together');
console.log('');

// Choose close primes
const p_close = 1009n; // Prime
const q_close = 1013n; // Next prime
const n_close = p_close * q_close;

console.log(`Close primes: p = ${p_close}, q = ${q_close}`);
console.log(`|p - q| = ${p_close > q_close ? p_close - q_close : q_close - p_close}`);
console.log(`n = p * q = ${n_close}`);

// Fermat factorization
function fermatFactor(n: bigint): [bigint, bigint] {
  let a = isqrt(n);
  if (a * a < n) a += 1n;

  let iterations = 0n;

  while (true) {
    const b2 = a * a - n;
    const b = isqrt(b2);

    iterations++;

    if (b * b === b2) {
      // Found it!
      const p = a + b;
      const q = a - b;
      console.log(`Fermat factorization succeeded in ${iterations} iterations`);
      return [p, q];
    }

    a += 1n;

    // Safety limit
    if (iterations > 10000n) {
      throw new Error('Too many iterations');
    }
  }
}

console.log('\nAttempting Fermat factorization...');
const [p_found, q_found] = fermatFactor(n_close);
console.log(`Found factors: p = ${p_found}, q = ${q_found}`);
console.log(`Verification: p * q = ${p_found * q_found}`);
console.log(`Correct: ${p_found * q_found === n_close}`);

// Compare with distant primes
console.log('\n---');
const p_distant = 2n ** 10n - 3n; // 1021, prime
const q_distant = 2n ** 9n + 9n; // 521, prime
const n_distant = p_distant * q_distant;

console.log(`\nDistant primes: p = ${p_distant}, q = ${q_distant}`);
console.log(`|p - q| = ${p_distant - q_distant}`);
console.log(`n = ${n_distant}`);

// This would take much longer
console.log('Fermat factorization would take many more iterations...');

console.log('\nMitigation:');
console.log('  - Ensure |p - q| > 2^(n_bits/2 - 100)');
console.log('  - Choose p and q independently and randomly');

// ============================================================================
// Attack 5: Low Private Exponent (Wiener's Attack)
// ============================================================================

console.log("\n=== Attack 5: Wiener's Attack (Low Private Exponent) ===\n");

/**
 * Wiener's attack applies when d < (1/3) * n^(1/4).
 * It uses continued fractions to find d from e and n.
 *
 * The idea: e*d = 1 + k*phi(n) for some k.
 * So e/phi(n) ≈ k/d.
 * Since phi(n) ≈ n, we have e/n ≈ k/d.
 * The convergents of the continued fraction of e/n include k/d.
 */

console.log('Scenario: Private exponent d is very small');
console.log('');

// This is a demonstration of the concept
// In practice, Wiener's attack requires careful continued fraction computation

const n_wiener = 3233n;
const phi_wiener = 3120n; // (61-1)(53-1)

// Choose a small d (vulnerable)
const d_small = 7n; // Very small d
// Compute corresponding e
const e_wiener = inverse_mod(d_small, phi_wiener);

console.log(`n = ${n_wiener}`);
console.log(`phi(n) = ${phi_wiener}`);
console.log(`Small d = ${d_small}`);
console.log(`Corresponding e = ${e_wiener}`);

// Check Wiener's bound
const n_fourth_root = isqrt(isqrt(n_wiener));
const wiener_bound = n_fourth_root / 3n;
console.log(`\nWiener's bound: d < n^(1/4)/3 ≈ ${wiener_bound}`);
console.log(`d = ${d_small} < ${wiener_bound}? ${d_small < wiener_bound || d_small <= 3n}`);

console.log('\nSimplified attack concept:');
console.log('  1. Compute continued fraction expansion of e/n');
console.log('  2. For each convergent k/d\':');
console.log('     a. Check if d\' gives valid phi(n)');
console.log('     b. If phi(n) = n - (n - phi(n) + 1) + 1 factors correctly, done');

// Demonstrate the relationship
console.log(`\nDemonstrating the relationship:`);
console.log(`  e/n = ${e_wiener}/${n_wiener}`);
console.log(`  e*d = ${e_wiener * d_small}`);
console.log(`  e*d mod phi(n) = ${(e_wiener * d_small) % phi_wiener}`);
console.log(`  k = (e*d - 1) / phi(n) = ${(e_wiener * d_small - 1n) / phi_wiener}`);

console.log('\nMitigation:');
console.log('  - Always use d > n^(1/4)');
console.log('  - Standard key generation (d ≈ n) is safe');
console.log('  - Use CRT representation instead of small d for speed');

// ============================================================================
// Attack 6: Timing Attack (Concept)
// ============================================================================

console.log('\n=== Attack 6: Timing Attack (Concept) ===\n');

/**
 * If RSA implementation time depends on the private key bits,
 * an attacker can recover d by measuring decryption times.
 *
 * Square-and-multiply has different timing based on whether
 * each bit of d is 0 or 1:
 * - Bit = 0: Only squaring
 * - Bit = 1: Squaring AND multiplication
 */

console.log('Timing attack concept:');
console.log('');
console.log('  Standard square-and-multiply for c^d:');
console.log('    result = 1');
console.log('    for each bit b of d (LSB to MSB):');
console.log('      if b == 1:');
console.log('        result = result * c  // EXTRA operation when bit is 1');
console.log('      c = c * c');
console.log('');
console.log('  Attack:');
console.log('    1. Send many ciphertexts');
console.log('    2. Measure decryption time for each');
console.log('    3. Statistical analysis reveals which bits of d are 1');

// Demonstrate timing difference
console.log('\nDemonstrating timing variation:');

const c_timing = 1234n;
const d_timing = 2753n;
const n_timing = 3233n;

// Simulate counting operations
function countOperations(c: bigint, d: bigint, n: bigint): number {
  let ops = 0;
  let dCopy = d;

  while (dCopy > 0n) {
    if (dCopy & 1n) {
      ops += 2; // Square + multiply
    } else {
      ops += 1; // Just square
    }
    dCopy >>= 1n;
  }

  return ops;
}

const d_test1 = 0b111111111n; // All 1s
const d_test2 = 0b100000000n; // One 1

console.log(`  d = ${d_test1} (binary: ${d_test1.toString(2)}):`);
console.log(`    Operations: ${countOperations(c_timing, d_test1, n_timing)}`);

console.log(`  d = ${d_test2} (binary: ${d_test2.toString(2)}):`);
console.log(`    Operations: ${countOperations(c_timing, d_test2, n_timing)}`);

console.log('\nMitigation:');
console.log('  - Constant-time implementation (always same number of operations)');
console.log('  - RSA blinding (randomize input)');
console.log('  - Montgomery ladder algorithm');

// ============================================================================
// Attack 7: Padding Oracle Attack (Bleichenbacher)
// ============================================================================

console.log('\n=== Attack 7: Padding Oracle Attack ===\n');

/**
 * If a server reveals whether PKCS#1 v1.5 padding is valid,
 * an attacker can decrypt arbitrary ciphertexts.
 *
 * PKCS#1 v1.5 padding: 0x00 || 0x02 || [random non-zero bytes] || 0x00 || [message]
 *
 * Attack:
 * 1. Multiply ciphertext by s^e mod n for various s
 * 2. Server reveals if decryption has valid padding
 * 3. Use this oracle to narrow down the plaintext range
 * 4. After ~1 million queries, plaintext is recovered
 */

console.log('Bleichenbacher\'s attack concept:');
console.log('');
console.log('  PKCS#1 v1.5 padding structure:');
console.log('    0x00 || 0x02 || [random bytes] || 0x00 || [message]');
console.log('');
console.log('  Oracle reveals: "Does decryption start with 0x00 0x02?"');
console.log('');
console.log('  Attack uses RSA malleability:');
console.log('    c\' = c * s^e mod n');
console.log('    Decrypts to: (m * s) mod n');
console.log('');
console.log('  Binary search using oracle responses:');
console.log('    - If padding valid: m*s is in "valid" range');
console.log('    - Narrow down m through ~2^20 queries');

// Demonstrate malleability exploitation
console.log('\nDemonstrating malleability:');

const m_pkcs = 100n;
const n_pkcs = 3233n;
const e_pkcs = 17n;
const c_pkcs = Mod(m_pkcs, n_pkcs).pow(e_pkcs).value;

console.log(`  Original: m = ${m_pkcs}, c = ${c_pkcs}`);

// Attacker modifies ciphertext
const s_attack = 2n;
const s_encrypted = Mod(s_attack, n_pkcs).pow(e_pkcs).value;
const c_modified = (c_pkcs * s_encrypted) % n_pkcs;

console.log(`  Modified: c' = c * ${s_attack}^e mod n = ${c_modified}`);

// If decrypted (by server), gives m * s
// Attacker learns information about m from padding oracle

console.log('\nMitigation:');
console.log('  - Use RSA-OAEP instead of PKCS#1 v1.5');
console.log('  - Never reveal padding validity directly');
console.log('  - Constant-time padding verification');

// ============================================================================
// Summary: Defense Checklist
// ============================================================================

console.log('\n=== RSA Security Checklist ===\n');

console.log('KEY GENERATION:');
console.log('  [ ] Use properly random prime generation');
console.log('  [ ] Ensure |p - q| is large (> 2^(n/2 - 100))');
console.log('  [ ] Check p-1 and q-1 have large prime factors');
console.log('  [ ] Use minimum 2048-bit keys');
console.log('');

console.log('PUBLIC EXPONENT:');
console.log('  [ ] Use e = 65537 (not e = 3)');
console.log('  [ ] Verify gcd(e, phi(n)) = 1');
console.log('');

console.log('PADDING:');
console.log('  [ ] Use RSA-OAEP for encryption');
console.log('  [ ] Use RSA-PSS for signatures');
console.log('  [ ] NEVER use textbook RSA');
console.log('');

console.log('IMPLEMENTATION:');
console.log('  [ ] Constant-time operations');
console.log('  [ ] RSA blinding');
console.log('  [ ] Secure memory handling');
console.log('');

console.log('KEY MANAGEMENT:');
console.log('  [ ] Unique n for each key pair');
console.log('  [ ] Separate keys for signing/encryption');
console.log('  [ ] Proper key storage (HSM when possible)');

console.log('\n' + '='.repeat(60));
console.log('RSA Attacks Tutorial Complete');
console.log('='.repeat(60));
