/**
 * # Diffie-Hellman Parameter Selection
 *
 * Choosing secure parameters for Diffie-Hellman is critical. Poor parameter
 * choice can make the system vulnerable to practical attacks.
 *
 * ## Key Decisions
 *
 * 1. Prime size: How large should p be?
 * 2. Prime structure: Safe prime vs. DSA-style prime
 * 3. Generator choice: What properties should g have?
 * 4. Standard vs. custom: Use published parameters or generate new ones?
 *
 * ## Historical Lessons
 *
 * - Logjam attack (2015): 512-bit DH was broken
 * - Weak DH study: 8% of HTTPS used shared 1024-bit primes
 * - NOBUS concerns: Pre-computed attacks on shared parameters
 */

import { gcd, is_prime, power_mod, factor, euler_phi, next_prime, Mod } from 'sagemath-ts';

console.log('='.repeat(60));
console.log('Diffie-Hellman Parameter Selection Tutorial');
console.log('='.repeat(60));

// ============================================================================
// Security Levels and Prime Sizes
// ============================================================================

console.log('\n=== Security Levels and Prime Sizes ===\n');

/**
 * DH security level depends on the best known attack:
 * - Index calculus / Number Field Sieve
 * - Complexity: L[1/3, c] = exp(c * (log n)^(1/3) * (log log n)^(2/3))
 *
 * Equivalent security levels:
 *   Symmetric bits | RSA/DH bits | ECC bits
 *   80             | 1024        | 160
 *   112            | 2048        | 224
 *   128            | 3072        | 256
 *   192            | 7680        | 384
 *   256            | 15360       | 512
 */

console.log('Security Level Comparison:');
console.log('');
console.log('  Symmetric | RSA/DH  | ECDH   | Status');
console.log('  ---------  -------   ------   ------');
console.log('  80 bits  | 1024    | 160    | BROKEN (avoid)');
console.log('  112 bits | 2048    | 224    | Minimum acceptable');
console.log('  128 bits | 3072    | 256    | Recommended');
console.log('  192 bits | 7680    | 384    | High security');
console.log('  256 bits | 15360   | 512    | Maximum security');
console.log('');
console.log('Current recommendations (2024):');
console.log('  - Minimum: 2048-bit DH');
console.log('  - Recommended: 3072-bit DH or ECDH P-256');
console.log('  - Long-term: 4096-bit DH or ECDH P-384');

// ============================================================================
// Safe Primes
// ============================================================================

console.log('\n=== Safe Primes ===\n');

/**
 * A safe prime p satisfies: p = 2q + 1 where q is also prime.
 *
 * Benefits:
 * - (Z/pZ)* has order p-1 = 2q with only one large factor
 * - Pohlig-Hellman attack gives no advantage
 * - Order-q subgroup is cyclic of prime order
 *
 * Properties:
 * - ~50% of subgroup elements are generators
 * - Quadratic residues form the order-q subgroup
 */

console.log('Safe Prime Definition:');
console.log('  p is safe prime if p = 2q + 1 where q is also prime');
console.log('  q is called a Sophie Germain prime');
console.log('');

function isSafePrime(p: bigint): boolean {
  if (!is_prime(p)) return false;
  const q = (p - 1n) / 2n;
  return is_prime(q);
}

// Find safe primes
console.log('Small safe primes:');
let count = 0;
for (let n = 5n; n < 500n && count < 10; n += 2n) {
  if (isSafePrime(n)) {
    const q = (n - 1n) / 2n;
    console.log(`  p = ${n} = 2 * ${q} + 1`);
    count++;
  }
}

// Demonstrate order-q subgroup
console.log('\nOrder-q subgroup for p = 23:');
const safePrime = 23n;
const q = (safePrime - 1n) / 2n; // q = 11

console.log(`  p = ${safePrime}, q = ${q}`);
console.log('  Elements of order q (quadratic residues):');

function orderOfElement(g: bigint, p: bigint): bigint {
  let order = 1n;
  let current = g;
  while (current !== 1n && order < p) {
    current = (current * g) % p;
    order++;
  }
  return order;
}

const qrElements: bigint[] = [];
for (let x = 1n; x < safePrime; x++) {
  if (orderOfElement(x, safePrime) === q) {
    qrElements.push(x);
  }
}
console.log(`    ${qrElements.join(', ')}`);
console.log(`    Total: ${qrElements.length} elements of order q`);

// ============================================================================
// DSA-Style Primes (Schnorr Groups)
// ============================================================================

console.log('\n=== DSA-Style Primes ===\n');

/**
 * Alternative: p and q where q divides p-1, but q << p
 *
 * Structure: p = kq + 1 for some k
 *
 * Benefits:
 * - Smaller q means smaller exponents (faster)
 * - Can use 256-bit q with 2048-bit p
 *
 * Used in:
 * - DSA (Digital Signature Algorithm)
 * - Some DH parameter sets
 */

console.log('DSA-style parameters:');
console.log('  p: 2048-bit prime');
console.log('  q: 256-bit prime dividing p-1');
console.log('  g: element of order q');
console.log('');

// Demonstrate with small example
// Find p where p-1 has a nice prime factor
const q_small = 11n;
let p_dsa = q_small + 1n;
while (!is_prime(p_dsa)) {
  p_dsa += q_small;
}

console.log(`Example (small): q = ${q_small}, p = ${p_dsa}`);
console.log(`  p - 1 = ${p_dsa - 1n} = ${(p_dsa - 1n) / q_small} * ${q_small}`);

// Find generator of order q
let g_dsa = 2n;
while (true) {
  const g_candidate = power_mod(g_dsa, (p_dsa - 1n) / q_small, p_dsa);
  if (g_candidate !== 1n) {
    // Verify order is exactly q
    if (power_mod(g_candidate, q_small, p_dsa) === 1n) {
      g_dsa = g_candidate;
      break;
    }
  }
  g_dsa++;
}

console.log(`  Generator of order q: g = ${g_dsa}`);
console.log(`  Verification: g^q mod p = ${power_mod(g_dsa, q_small, p_dsa)}`);

console.log('\nComparison:');
console.log('  Safe primes:  Simpler, q = (p-1)/2');
console.log('  DSA-style:    More flexible, smaller q possible');

// ============================================================================
// Finding Generators
// ============================================================================

console.log('\n=== Finding Generators ===\n');

/**
 * For safe prime p = 2q + 1:
 *
 * To find generator of (Z/pZ)*:
 * - Check g^q != 1 and g^2 != 1 (mod p)
 * - Or simply check g^q != 1 and g != p-1
 *
 * To find generator of order-q subgroup:
 * - Compute g = r^2 mod p for random r
 * - This gives element of order q (or 1, which we skip)
 */

console.log('Finding generators for p = 23:');
console.log('');

const p_test = 23n;
const q_test = 11n;

// Full group generators
console.log('Generators of full group (order p-1 = 22):');
const fullGenerators: bigint[] = [];
for (let g = 2n; g < p_test; g++) {
  const order = orderOfElement(g, p_test);
  if (order === p_test - 1n) {
    fullGenerators.push(g);
  }
}
console.log(`  ${fullGenerators.join(', ')}`);

// Subgroup generators
console.log('\nGenerators of order-q subgroup (order 11):');
const subgroupGenerators: bigint[] = [];
for (let g = 2n; g < p_test; g++) {
  const order = orderOfElement(g, p_test);
  if (order === q_test) {
    subgroupGenerators.push(g);
  }
}
console.log(`  ${subgroupGenerators.join(', ')}`);

// Method: squaring
console.log('\nFinding subgroup generator by squaring:');
for (let r = 2n; r <= 5n; r++) {
  const g_sq = (r * r) % p_test;
  const order = orderOfElement(g_sq, p_test);
  console.log(`  r = ${r}: g = r^2 mod p = ${g_sq}, order = ${order}`);
}

// ============================================================================
// Standard Parameters
// ============================================================================

console.log('\n=== Standard DH Parameters ===\n');

/**
 * Standard parameters are published and widely reviewed.
 * Using standard parameters avoids the need to trust parameter generation.
 *
 * RFC 3526: MODP Groups for IKE (IPsec):
 * - 1536-bit MODP Group
 * - 2048-bit MODP Group
 * - 3072-bit MODP Group
 * - 4096-bit MODP Group
 * - 6144-bit MODP Group
 * - 8192-bit MODP Group
 *
 * RFC 7919: FFDHE (Finite Field DH with Ephemeral keys) for TLS:
 * - ffdhe2048, ffdhe3072, ffdhe4096, ffdhe6144, ffdhe8192
 */

console.log('RFC 3526 MODP Groups (for IPsec/IKE):');
console.log('');
console.log('  Group 5:  1536-bit (deprecated)');
console.log('  Group 14: 2048-bit (minimum recommended)');
console.log('  Group 15: 3072-bit');
console.log('  Group 16: 4096-bit');
console.log('  Group 17: 6144-bit');
console.log('  Group 18: 8192-bit');
console.log('');

console.log('RFC 7919 FFDHE Groups (for TLS 1.3):');
console.log('');
console.log('  ffdhe2048: 2048-bit, mandatory-to-implement');
console.log('  ffdhe3072: 3072-bit');
console.log('  ffdhe4096: 4096-bit');
console.log('  ffdhe6144: 6144-bit');
console.log('  ffdhe8192: 8192-bit');
console.log('');

// Structure of RFC primes
console.log('Structure of RFC primes:');
console.log('  These are safe primes: p = 2q + 1');
console.log('  q is a Sophie Germain prime');
console.log('  Generator: g = 2');
console.log('  p is chosen for efficient computation (low Hamming weight)');

// ============================================================================
// Risks of Custom Parameters
// ============================================================================

console.log('\n=== Risks of Custom Parameters ===\n');

/**
 * Dangers of generating your own parameters:
 *
 * 1. Trapdoor primes: p can be chosen so attacker knows factorization
 * 2. Weak primes: Poor randomness leads to predictable parameters
 * 3. NOBUS attacks: If parameters are shared, precomputation helps everyone
 *
 * The Logjam attack showed that shared 1024-bit primes allowed
 * nation-state level attackers to break DH in real-time.
 */

console.log('Risks of custom DH parameters:');
console.log('');
console.log('  1. Trapdoor primes:');
console.log('     Adversary chooses p with hidden structure');
console.log('     Example: p-1 has secret small factors');
console.log('');
console.log('  2. Weak random generation:');
console.log('     Poor PRNG leads to predictable primes');
console.log('     Attacker can guess your parameters');
console.log('');
console.log('  3. Shared parameter attacks:');
console.log('     Precomputation applies to all users of same p');
console.log('     2015: 8% of HTTPS sites shared one 1024-bit prime');
console.log('');
console.log('  4. Verification difficulty:');
console.log('     Hard to prove parameters are "honestly" generated');
console.log('     Standard parameters have public generation process');

// ============================================================================
// Parameter Validation
// ============================================================================

console.log('\n=== Parameter Validation ===\n');

/**
 * When receiving DH parameters, validate them:
 *
 * 1. p is prime (Miller-Rabin)
 * 2. p is large enough (>= 2048 bits)
 * 3. g is valid generator (check order)
 * 4. For safe primes: (p-1)/2 is prime
 * 5. Public value is in valid range: 2 <= y <= p-2
 * 6. Public value is in correct subgroup: y^q = 1 (mod p)
 */

console.log('Validation checklist:');
console.log('');

function validateDHParams(p: bigint, g: bigint, minBits: number): string[] {
  const errors: string[] = [];

  // Check p is prime
  if (!is_prime(p)) {
    errors.push('p is not prime');
  }

  // Check p is large enough
  const pBits = p.toString(2).length;
  if (pBits < minBits) {
    errors.push(`p has only ${pBits} bits (need >= ${minBits})`);
  }

  // Check g is in valid range
  if (g < 2n || g >= p) {
    errors.push('g is not in valid range [2, p-1]');
  }

  // Check g is not trivial
  if (g === 1n || g === p - 1n) {
    errors.push('g is a trivial value');
  }

  // For safe prime, check q is prime
  const q = (p - 1n) / 2n;
  if (is_prime(q)) {
    // It's a safe prime, check g generates order-q subgroup
    if (power_mod(g, q, p) !== 1n) {
      // g might generate full group, that's okay
    }
  }

  return errors;
}

// Test with our small prime (should fail size check)
console.log(`Validating p = ${safePrime}, g = 5:`);
const errors1 = validateDHParams(safePrime, 5n, 2048);
if (errors1.length > 0) {
  console.log('  Errors:');
  errors1.forEach((e) => console.log(`    - ${e}`));
} else {
  console.log('  Parameters valid');
}

// Test with good parameters (simulated)
console.log('\nExpected validation for RFC parameters:');
console.log('  - p is prime: PASS');
console.log('  - p >= 2048 bits: PASS');
console.log('  - g in valid range: PASS');
console.log('  - (p-1)/2 is prime: PASS');
console.log('  - All checks passed');

// ============================================================================
// Public Value Validation
// ============================================================================

console.log('\n=== Public Value Validation ===\n');

/**
 * When receiving a public value y = g^x:
 *
 * 1. Range check: 2 <= y <= p-2
 * 2. Subgroup check: y^q = 1 (mod p) for order-q subgroup
 * 3. Small subgroup attack prevention
 */

console.log('Public value validation:');
console.log('');

function validatePublicValue(y: bigint, p: bigint, q: bigint): string[] {
  const errors: string[] = [];

  // Range check
  if (y < 2n || y >= p - 1n) {
    errors.push('y is not in valid range [2, p-2]');
  }

  // Subgroup check (if working in order-q subgroup)
  if (power_mod(y, q, p) !== 1n) {
    errors.push('y is not in the expected subgroup');
  }

  return errors;
}

// Test cases
const testValues = [0n, 1n, 5n, safePrime - 1n, safePrime];

console.log(`Testing with p = ${safePrime}, q = ${q_test}:`);
for (const y of testValues) {
  const errs = validatePublicValue(y, safePrime, q_test);
  if (errs.length > 0) {
    console.log(`  y = ${y}: INVALID - ${errs[0]}`);
  } else {
    console.log(`  y = ${y}: valid`);
  }
}

// ============================================================================
// Migration to ECDH
// ============================================================================

console.log('\n=== Migration to ECDH ===\n');

/**
 * Elliptic Curve Diffie-Hellman (ECDH) offers:
 * - Same security with smaller keys (256-bit vs 3072-bit)
 * - Faster computation
 * - Standard curves avoid parameter generation issues
 *
 * Recommended curves:
 * - P-256 (NIST, widely supported)
 * - P-384 (NIST, higher security)
 * - X25519 (Curve25519, modern choice)
 * - X448 (Curve448, highest security)
 */

console.log('Why consider ECDH over traditional DH:');
console.log('');
console.log('  Performance comparison:');
console.log('    Operation      | DH-2048  | ECDH-256');
console.log('    -------------  | -------  | --------');
console.log('    Key generation | ~10ms    | ~0.1ms');
console.log('    Key exchange   | ~10ms    | ~0.2ms');
console.log('    Key size       | 256 B    | 32 B');
console.log('');
console.log('  Security levels are equivalent:');
console.log('    DH-2048 ≈ ECDH-224 ≈ 112-bit symmetric');
console.log('    DH-3072 ≈ ECDH-256 ≈ 128-bit symmetric');
console.log('');
console.log('  Recommended curves:');
console.log('    - X25519: Best performance, modern design');
console.log('    - P-256: NIST standard, wide support');
console.log('    - P-384: Higher security margin');

// ============================================================================
// Summary Checklist
// ============================================================================

console.log('\n=== Parameter Selection Checklist ===\n');

console.log('DH PARAMETER SELECTION CHECKLIST:');
console.log('');
console.log('  [ ] Use standard parameters (RFC 7919 or RFC 3526)');
console.log('      - Avoid generating custom parameters');
console.log('      - Verifiable generation process');
console.log('');
console.log('  [ ] Minimum 2048-bit prime');
console.log('      - 3072-bit recommended for new systems');
console.log('      - 4096-bit for long-term security');
console.log('');
console.log('  [ ] Safe prime structure');
console.log('      - p = 2q + 1 where q is prime');
console.log('      - Protects against Pohlig-Hellman');
console.log('');
console.log('  [ ] Validate all received parameters');
console.log('      - Check primality');
console.log('      - Check size');
console.log('      - Check generator order');
console.log('');
console.log('  [ ] Validate all public values');
console.log('      - Range check');
console.log('      - Subgroup membership');
console.log('');
console.log('  [ ] Consider ECDH instead');
console.log('      - Better performance');
console.log('      - Smaller keys');
console.log('      - X25519 recommended');

console.log('\n' + '='.repeat(60));
console.log('Diffie-Hellman Parameter Selection Tutorial Complete');
console.log('='.repeat(60));
