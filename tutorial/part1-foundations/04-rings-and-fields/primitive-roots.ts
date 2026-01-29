/**
 * # Primitive Roots: Generators of (Z/pZ)*
 *
 * A primitive root modulo n is a generator of the multiplicative group (Z/nZ)*.
 * For prime p, primitive roots are essential for Diffie-Hellman key exchange
 * and other cryptographic protocols.
 *
 * ## Mathematical Background
 *
 * An element g in (Z/nZ)* is a PRIMITIVE ROOT if ord(g) = phi(n).
 * Equivalently, g generates the entire group: (Z/nZ)* = <g>.
 *
 * Key facts:
 *   - Primitive roots exist iff n = 1, 2, 4, p^k, or 2p^k (p odd prime)
 *   - For prime p, there are exactly phi(p-1) primitive roots
 *   - If g is a primitive root, then g^k is a primitive root iff gcd(k, phi(n)) = 1
 *
 * ## Why This Matters for Cryptography
 *
 *   - Diffie-Hellman: The generator g must be a primitive root (or generate a large subgroup)
 *   - Discrete Log Problem: Hardness depends on working in a large cyclic group
 *   - ElGamal encryption: Requires a primitive root
 *   - Schnorr signatures: Use a prime-order subgroup
 */

import { gcd, euler_phi, power_mod, factor, divisors, is_prime, inverse_mod, Zmod } from 'sagemath-ts';

// =============================================================================
// Example 1: Finding Primitive Roots by Brute Force
// =============================================================================
console.log("=== Example 1: Finding Primitive Roots ===\n");

/**
 * A primitive root g has order phi(n).
 * We can find one by checking each candidate.
 */
function computeOrder(g: bigint, n: bigint): bigint {
  if (gcd(g, n) !== 1n) return 0n;

  let order = 1n;
  let current = g % n;
  const maxOrder = euler_phi(n);

  while (current !== 1n && order <= maxOrder) {
    current = (current * g) % n;
    order++;
  }

  return current === 1n ? order : 0n;
}

function findPrimitiveRoots(n: bigint): bigint[] {
  const phi = euler_phi(n);
  const roots: bigint[] = [];

  for (let g = 2n; g < n; g++) {
    if (gcd(g, n) === 1n && computeOrder(g, n) === phi) {
      roots.push(g);
    }
  }

  return roots;
}

const p = 11n;
const primitiveRoots = findPrimitiveRoots(p);

console.log(`Primitive roots modulo ${p}:`);
console.log(`  (Z/${p}Z)* has order ${p - 1n}`);
console.log(`  Primitive roots: {${primitiveRoots.join(', ')}}`);
console.log(`  Count: ${primitiveRoots.length} = phi(${p - 1n}) = ${euler_phi(p - 1n)}`);

console.log("\nVerifying orders:");
for (const g of primitiveRoots) {
  console.log(`  ord(${g}) = ${computeOrder(g, p)}`);
}

console.log();

// =============================================================================
// Example 2: Efficient Primitive Root Test
// =============================================================================
console.log("=== Example 2: Efficient Primitive Root Test ===\n");

/**
 * Instead of computing the full order, we can test if g is a primitive root
 * by checking that g^((p-1)/q) != 1 for each prime factor q of p-1.
 */
function isPrimitiveRoot(g: bigint, p: bigint): boolean {
  if (!is_prime(p) || gcd(g, p) !== 1n) return false;

  const pMinus1 = p - 1n;
  const primeFactors = factor(pMinus1)
    .filter(([f]) => f > 0n)
    .map(([f]) => f);

  for (const q of primeFactors) {
    if (power_mod(g, pMinus1 / q, p) === 1n) {
      return false;
    }
  }

  return true;
}

const testP = 23n;
const pMinus1 = testP - 1n;
const primeFactors = factor(pMinus1).filter(([f]) => f > 0n).map(([f]) => f);

console.log(`Testing primitive roots modulo ${testP}:`);
console.log(`  ${testP} - 1 = ${pMinus1} = ${factor(pMinus1).filter(([f]) => f > 0n).map(([p, e]) => e === 1n ? `${p}` : `${p}^${e}`).join(' * ')}`);
console.log(`  Prime factors to check: {${primeFactors.join(', ')}}\n`);

for (let g = 2n; g <= 10n; g++) {
  const isPrim = isPrimitiveRoot(g, testP);
  console.log(`  g = ${g}:`);

  for (const q of primeFactors) {
    const exp = pMinus1 / q;
    const result = power_mod(g, exp, testP);
    console.log(`    g^(${pMinus1}/${q}) = g^${exp} = ${result} ${result === 1n ? '= 1 (FAIL)' : '!= 1 (ok)'}`);
  }

  console.log(`    Primitive root? ${isPrim}\n`);
}

// =============================================================================
// Example 3: All Primitive Roots from One
// =============================================================================
console.log("=== Example 3: Generating All Primitive Roots ===\n");

/**
 * If g is a primitive root mod p, then g^k is also a primitive root
 * iff gcd(k, p-1) = 1.
 */
const genP = 13n;
const genOrder = genP - 1n;

// Find smallest primitive root
let firstRoot = 2n;
while (!isPrimitiveRoot(firstRoot, genP)) firstRoot++;

console.log(`For p = ${genP}:`);
console.log(`  First primitive root: g = ${firstRoot}`);
console.log(`  Order = ${genOrder}`);
console.log(`\nAll primitive roots are g^k where gcd(k, ${genOrder}) = 1:\n`);

const allRoots: bigint[] = [];
for (let k = 1n; k < genOrder; k++) {
  if (gcd(k, genOrder) === 1n) {
    const root = power_mod(firstRoot, k, genP);
    allRoots.push(root);
    console.log(`  g^${k} = ${firstRoot}^${k} = ${root} mod ${genP}`);
  }
}

console.log(`\nAll primitive roots: {${allRoots.sort((a, b) => Number(a - b)).join(', ')}}`);
console.log(`Count: ${allRoots.length} = phi(${genOrder}) = ${euler_phi(genOrder)}`);

console.log();

// =============================================================================
// Example 4: Discrete Logarithm
// =============================================================================
console.log("=== Example 4: Discrete Logarithm ===\n");

/**
 * If g is a primitive root mod p, every element a in (Z/pZ)* can be written
 * as a = g^x for some x. The value x is the discrete logarithm of a base g.
 */
const dlogP = 11n;
const dlogG = findPrimitiveRoots(dlogP)[0]!;

console.log(`Discrete logarithms in (Z/${dlogP}Z)* with generator g = ${dlogG}:\n`);

// Build discrete log table
const dlogTable = new Map<bigint, bigint>();
let power = 1n;
for (let x = 0n; x < dlogP - 1n; x++) {
  dlogTable.set(power, x);
  power = (power * dlogG) % dlogP;
}

console.log("Discrete log table:");
for (let a = 1n; a < dlogP; a++) {
  const x = dlogTable.get(a)!;
  console.log(`  log_${dlogG}(${a}) = ${x}  (verify: ${dlogG}^${x} = ${power_mod(dlogG, x, dlogP)})`);
}

console.log("\nThe discrete log problem: Given g, p, and h = g^x, find x.");
console.log("This is computationally hard for large primes!");

console.log();

// =============================================================================
// Example 5: Index Calculus Preview
// =============================================================================
console.log("=== Example 5: Discrete Log Properties ===\n");

/**
 * Discrete logarithms satisfy:
 *   log_g(ab) = log_g(a) + log_g(b) (mod p-1)
 *   log_g(a^k) = k * log_g(a) (mod p-1)
 */
const icP = 13n;
const icG = findPrimitiveRoots(icP)[0]!;
const icOrder = icP - 1n;

// Build dlog table
const icDlog = new Map<bigint, bigint>();
let icPower = 1n;
for (let x = 0n; x < icOrder; x++) {
  icDlog.set(icPower, x);
  icPower = (icPower * icG) % icP;
}

console.log(`Discrete log properties in (Z/${icP}Z)* with g = ${icG}:\n`);

// Test multiplicative property
const a = 3n, b = 5n;
const ab = (a * b) % icP;
const logA = icDlog.get(a)!, logB = icDlog.get(b)!, logAB = icDlog.get(ab)!;

console.log("Multiplicative property: log(ab) = log(a) + log(b) mod (p-1)");
console.log(`  a = ${a}, b = ${b}, ab = ${ab}`);
console.log(`  log(${a}) = ${logA}`);
console.log(`  log(${b}) = ${logB}`);
console.log(`  log(${a}) + log(${b}) = ${logA + logB} = ${(logA + logB) % icOrder} mod ${icOrder}`);
console.log(`  log(${ab}) = ${logAB}`);
console.log(`  Equal? ${(logA + logB) % icOrder === logAB}\n`);

// Test power property
const k = 3n;
const ak = power_mod(a, k, icP);
const logAK = icDlog.get(ak)!;

console.log("Power property: log(a^k) = k * log(a) mod (p-1)");
console.log(`  a = ${a}, k = ${k}, a^k = ${ak}`);
console.log(`  k * log(${a}) = ${k} * ${logA} = ${(k * logA) % icOrder} mod ${icOrder}`);
console.log(`  log(${ak}) = ${logAK}`);
console.log(`  Equal? ${(k * logA) % icOrder === logAK}`);

console.log();

// =============================================================================
// Example 6: Primitive Roots for Composite Moduli
// =============================================================================
console.log("=== Example 6: Primitive Roots for Composite Moduli ===\n");

/**
 * Primitive roots exist for n = 2, 4, p^k, 2p^k (p odd prime).
 * Let's verify this for small moduli.
 */
function hasPrimitiveRoot(n: bigint): boolean {
  if (n === 1n || n === 2n || n === 4n) return true;

  const factors = factor(n).filter(([p]) => p > 0n);

  // n = p^k for odd prime p
  if (factors.length === 1) {
    const [prime, _exp] = factors[0]!;
    return prime !== 2n;
  }

  // n = 2 * p^k for odd prime p
  if (factors.length === 2) {
    const [[p1, e1], [p2, _e2]] = factors;
    if (p1 === 2n && e1 === 1n && p2! > 2n) return true;
  }

  return false;
}

console.log("Which moduli have primitive roots?\n");

for (let n = 2n; n <= 25n; n++) {
  const has = hasPrimitiveRoot(n);
  const phi = euler_phi(n);
  const roots = findPrimitiveRoots(n);

  if (has) {
    console.log(`  n=${n.toString().padStart(2)}: YES, phi=${phi.toString().padStart(2)}, primitive roots: {${roots.slice(0, 5).join(', ')}${roots.length > 5 ? '...' : ''}}`);
  } else {
    console.log(`  n=${n.toString().padStart(2)}: NO (group is not cyclic)`);
  }
}

console.log();

// =============================================================================
// Example 7: Lifting Primitive Roots
// =============================================================================
console.log("=== Example 7: Primitive Roots Modulo Prime Powers ===\n");

/**
 * If g is a primitive root mod p, then either g or g + p is a
 * primitive root mod p^2, and can be lifted to any p^k.
 */
const liftP = 7n;
const basePrimRoot = findPrimitiveRoots(liftP)[0]!;

console.log(`Starting with primitive root ${basePrimRoot} mod ${liftP}:\n`);

// Check if g is a primitive root mod p^2
function isPrimitiveRootPk(g: bigint, p: bigint, k: bigint): boolean {
  const pk = p ** k;
  const phi = pk - pk / p;  // phi(p^k) = p^(k-1) * (p-1)

  // Check g^(phi/q) != 1 for each prime factor q of phi
  const phiFactors = factor(phi).filter(([f]) => f > 0n).map(([f]) => f);

  for (const q of phiFactors) {
    if (power_mod(g, phi / q, pk) === 1n) {
      return false;
    }
  }

  return gcd(g, pk) === 1n;
}

for (let k = 1n; k <= 3n; k++) {
  const pk = liftP ** k;
  const phi = pk - pk / liftP;

  console.log(`Mod ${liftP}^${k} = ${pk} (phi = ${phi}):`);

  // Try g and g + p
  for (const candidate of [basePrimRoot, basePrimRoot + liftP]) {
    const isPrim = isPrimitiveRootPk(candidate, liftP, k);
    if (isPrim) {
      console.log(`  ${candidate} is a primitive root mod ${pk}`);
      break;
    }
  }
}

console.log();

// =============================================================================
// Example 8: Safe Primes and Generators
// =============================================================================
console.log("=== Example 8: Safe Primes and Generators ===\n");

/**
 * A safe prime is p = 2q + 1 where q is also prime.
 * For safe primes, (Z/pZ)* has only subgroups of order 1, 2, q, and 2q.
 * Any element != +/-1 generates either the full group or the subgroup of order q.
 */
function isSafePrime(p: bigint): boolean {
  if (!is_prime(p)) return false;
  const q = (p - 1n) / 2n;
  return is_prime(q);
}

// Find some safe primes
const safePrimes: bigint[] = [];
for (let p = 5n; safePrimes.length < 5; p += 2n) {
  if (isSafePrime(p)) {
    safePrimes.push(p);
  }
}

console.log(`Safe primes: {${safePrimes.join(', ')}}\n`);

// Analyze the structure
const safeP = safePrimes[2]!;  // 23
const safeQ = (safeP - 1n) / 2n;

console.log(`For safe prime p = ${safeP} = 2 * ${safeQ} + 1:`);
console.log(`  (Z/${safeP}Z)* has order ${safeP - 1n}`);
console.log(`  Subgroups have orders: 1, 2, ${safeQ}, ${safeP - 1n}\n`);

console.log("Element classification:");
for (let g = 2n; g < safeP; g++) {
  const ord = computeOrder(g, safeP);
  let type: string;
  if (ord === 1n) type = "identity";
  else if (ord === 2n) type = "order 2 (= -1)";
  else if (ord === safeQ) type = `order q (subgroup)`;
  else type = "generator";

  console.log(`  ${g.toString().padStart(2)}: order ${ord.toString().padStart(2)}, ${type}`);
}

console.log("\nFor Diffie-Hellman with safe primes:");
console.log("  - Use g that generates full group (order 2q)");
console.log("  - Or use g^2 to generate subgroup of order q (recommended)");

console.log();

// =============================================================================
// Example 9: Pohlig-Hellman Attack Preview
// =============================================================================
console.log("=== Example 9: Why Group Order Matters ===\n");

/**
 * The Pohlig-Hellman algorithm can solve the discrete log problem
 * efficiently if p-1 has only small prime factors.
 */
console.log("Pohlig-Hellman attack vulnerability:\n");

// A "bad" prime where p-1 is smooth
const badP = 97n;  // 97 - 1 = 96 = 2^5 * 3
const badFactors = factor(badP - 1n);

console.log(`"Bad" prime: p = ${badP}`);
console.log(`  p - 1 = ${badP - 1n} = ${badFactors.filter(([f]) => f > 0n).map(([p, e]) => e === 1n ? `${p}` : `${p}^${e}`).join(' * ')}`);
console.log(`  Largest prime factor: ${Math.max(...badFactors.filter(([f]) => f > 0n).map(([f]) => Number(f)))}`);
console.log(`  Vulnerable to Pohlig-Hellman attack!\n`);

// A "good" prime (safe prime)
const goodP = 23n;
const goodFactors = factor(goodP - 1n);

console.log(`"Good" prime (safe): p = ${goodP}`);
console.log(`  p - 1 = ${goodP - 1n} = ${goodFactors.filter(([f]) => f > 0n).map(([p, e]) => e === 1n ? `${p}` : `${p}^${e}`).join(' * ')}`);
console.log(`  Largest prime factor: ${Math.max(...goodFactors.filter(([f]) => f > 0n).map(([f]) => Number(f)))}`);
console.log(`  Resistant to Pohlig-Hellman (use q-subgroup).\n`);

console.log("For cryptographic security:");
console.log("  - Use safe primes p = 2q + 1");
console.log("  - Or use primes where p-1 = 2q for large prime q");
console.log("  - Modern protocols often use prime-order subgroups");

console.log();

// =============================================================================
// Exercise: Find Primitive Root of a Large Prime
// =============================================================================
console.log("=== Exercise: Find Primitive Root ===\n");

/**
 * EXERCISE: Find a primitive root modulo 41.
 */
const exerciseP = 41n;
const exerciseOrder = exerciseP - 1n;  // = 40

console.log(`Finding a primitive root modulo ${exerciseP}:`);
console.log(`  ${exerciseP} - 1 = ${exerciseOrder} = ${factor(exerciseOrder).filter(([f]) => f > 0n).map(([p, e]) => e === 1n ? `${p}` : `${p}^${e}`).join(' * ')}\n`);

// Find first primitive root
let exerciseRoot = 2n;
while (!isPrimitiveRoot(exerciseRoot, exerciseP)) {
  exerciseRoot++;
}

console.log(`First primitive root: ${exerciseRoot}`);

// Verify
const expectedOrder = exerciseOrder;
const actualOrder = computeOrder(exerciseRoot, exerciseP);
console.log(`  Expected order: ${expectedOrder}`);
console.log(`  Actual order: ${actualOrder}`);
console.log(`  Is primitive root: ${expectedOrder === actualOrder}`);

// Count all primitive roots
const allPrimRoots = findPrimitiveRoots(exerciseP);
console.log(`\nTotal primitive roots: ${allPrimRoots.length} = phi(${exerciseOrder}) = ${euler_phi(exerciseOrder)}`);

console.log();

// =============================================================================
// Cryptographic Application: Diffie-Hellman Setup
// =============================================================================
console.log("=== Cryptographic Application: Diffie-Hellman Generator Selection ===\n");

console.log("Diffie-Hellman requires careful generator selection:\n");

const dhP = 23n;  // Safe prime
const dhQ = (dhP - 1n) / 2n;  // = 11

console.log(`1. Choose a safe prime p = ${dhP} = 2 * ${dhQ} + 1\n`);

// Find a generator of the order-q subgroup
let dhG = 2n;
while (power_mod(dhG, dhQ, dhP) === 1n || power_mod(dhG, 2n, dhP) === 1n) {
  dhG++;
}
// Use g^2 to ensure we're in the q-subgroup
const dhGenQ = power_mod(dhG, 2n, dhP);

console.log(`2. Find generator of order-q subgroup:`);
console.log(`   g = ${dhG}^2 = ${dhGenQ} (order ${dhQ})\n`);

// Simulate key exchange
const alicePrivate = 7n;
const bobPrivate = 9n;

const alicePublic = power_mod(dhGenQ, alicePrivate, dhP);
const bobPublic = power_mod(dhGenQ, bobPrivate, dhP);

const sharedAlice = power_mod(bobPublic, alicePrivate, dhP);
const sharedBob = power_mod(alicePublic, bobPrivate, dhP);

console.log(`3. Key exchange:`);
console.log(`   Alice: a = ${alicePrivate}, A = g^a = ${alicePublic}`);
console.log(`   Bob:   b = ${bobPrivate}, B = g^b = ${bobPublic}`);
console.log(`   Shared: g^(ab) = ${sharedAlice}\n`);

console.log("Security comes from the difficulty of computing discrete logs in <g>!");
