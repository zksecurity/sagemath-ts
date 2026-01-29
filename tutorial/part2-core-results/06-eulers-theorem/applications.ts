/**
 * # Applications of Euler's Theorem
 *
 * Euler's theorem and the totient function have numerous applications
 * in cryptography and number theory:
 *
 * 1. RSA cryptosystem
 * 2. Computing multiplicative order
 * 3. Primitive roots
 * 4. Discrete logarithm
 * 5. Modular arithmetic optimizations
 *
 * @module tutorial/part2/euler-applications
 */

import {
  euler_phi,
  power_mod,
  gcd,
  factor,
  is_prime,
  inverse_mod,
  Zmod,
  formatFactorization,
  divisors,
} from 'sagemath-ts';

console.log("=".repeat(70));
console.log("APPLICATIONS OF EULER'S THEOREM");
console.log("=".repeat(70));

// ============================================================================
// Application 1: RSA Cryptosystem
// ============================================================================
console.log("\n=== Application 1: RSA Cryptosystem ===\n");

console.log("RSA relies on Euler's theorem: a^phi(n) = 1 (mod n) for gcd(a,n) = 1\n");

// Key generation
const p = 61n;
const q = 53n;
const n = p * q;
const phiN = (p - 1n) * (q - 1n); // = euler_phi(n) for n = pq
const e = 17n; // Public exponent (must be coprime to phi(n))
const d = inverse_mod(e, phiN); // Private exponent

console.log("Key Generation:");
console.log(`  p = ${p}, q = ${q} (secret primes)`);
console.log(`  n = p * q = ${n} (public modulus)`);
console.log(`  phi(n) = (p-1)(q-1) = ${phiN} (secret)`);
console.log(`  e = ${e} (public exponent, gcd(e, phi(n)) = ${gcd(e, phiN)})`);
console.log(`  d = e^(-1) mod phi(n) = ${d} (private exponent)`);
console.log(`  Verify: e * d mod phi(n) = ${(e * d) % phiN}`);

console.log(`\n  Public key: (n=${n}, e=${e})`);
console.log(`  Private key: (n=${n}, d=${d})`);

// Encryption and decryption
const message = 42n;
const ciphertext = power_mod(message, e, n);
const decrypted = power_mod(ciphertext, d, n);

console.log(`\nEncryption/Decryption:`);
console.log(`  Message: m = ${message}`);
console.log(`  Encrypt: c = m^e mod n = ${message}^${e} mod ${n} = ${ciphertext}`);
console.log(`  Decrypt: m = c^d mod n = ${ciphertext}^${d} mod ${n} = ${decrypted}`);

console.log(`\nWhy decryption works (Euler's theorem):`);
console.log(`  c^d = (m^e)^d = m^(ed) mod n`);
console.log(`  Since ed = 1 + k*phi(n) for some k:`);
console.log(`    m^(ed) = m * (m^phi(n))^k`);
console.log(`  By Euler's theorem, m^phi(n) = 1 mod n (when gcd(m,n) = 1)`);
console.log(`  So m^(ed) = m * 1^k = m mod n`);

// ============================================================================
// Application 2: Finding Multiplicative Order
// ============================================================================
console.log("\n=== Application 2: Finding Multiplicative Order ===\n");

console.log("The multiplicative order ord_n(a) is the smallest positive k");
console.log("such that a^k = 1 (mod n). By Euler's theorem, ord_n(a) | phi(n).\n");

function multiplicativeOrder(a: bigint, n: bigint): bigint | null {
  if (gcd(a, n) !== 1n) return null;

  const phi = euler_phi(n);
  const divs = divisors(phi);

  // Find smallest divisor d of phi such that a^d = 1
  for (const d of divs) {
    if (power_mod(a, d, n) === 1n) {
      return d;
    }
  }

  return phi; // Should always be reached
}

const nOrder = 13n;
const phiOrder = euler_phi(nOrder);
console.log(`In Z/${nOrder}Z, phi(${nOrder}) = ${phiOrder}`);
console.log(`Divisors of phi: {${divisors(phiOrder).join(", ")}}\n`);

console.log("Element | Order | Divides phi?");
console.log("--------+-------+-------------");
for (let a = 1n; a < nOrder; a++) {
  const order = multiplicativeOrder(a, nOrder);
  const divides = order !== null && phiOrder % order === 0n;
  console.log(`   ${a.toString().padStart(2)}   |  ${order?.toString().padStart(2) ?? "N/A"}   | ${divides}`);
}

// ============================================================================
// Application 3: Primitive Roots
// ============================================================================
console.log("\n=== Application 3: Primitive Roots ===\n");

console.log("A primitive root mod n is an element of order phi(n).");
console.log("It generates the entire multiplicative group (Z/nZ)*.\n");

function isPrimitiveRoot(g: bigint, n: bigint): boolean {
  const phi = euler_phi(n);
  if (gcd(g, n) !== 1n) return false;

  // g is primitive root iff g^(phi/p) != 1 for all prime divisors p of phi
  const primeFactors = factor(phi).filter(([p, _]) => p > 0n).map(([p, _]) => p);

  for (const p of primeFactors) {
    if (power_mod(g, phi / p, n) === 1n) {
      return false;
    }
  }

  return true;
}

function findPrimitiveRoot(n: bigint): bigint | null {
  for (let g = 2n; g < n; g++) {
    if (isPrimitiveRoot(g, n)) {
      return g;
    }
  }
  return null;
}

const primesToTest = [7n, 11n, 13n, 23n];

for (const p of primesToTest) {
  const primitiveRoots: bigint[] = [];
  for (let g = 2n; g < p; g++) {
    if (isPrimitiveRoot(g, p)) {
      primitiveRoots.push(g);
    }
  }

  console.log(`p = ${p.toString().padStart(2)}: phi(${p}) = ${euler_phi(p)}`);
  console.log(`  Primitive roots: {${primitiveRoots.join(", ")}}`);
  console.log(`  Number of primitive roots: ${primitiveRoots.length} = phi(phi(${p})) = phi(${euler_phi(p)}) = ${euler_phi(euler_phi(p))}\n`);
}

console.log("Theorem: The number of primitive roots mod p equals phi(p-1).");

// ============================================================================
// Application 4: Order Finding for Shor's Algorithm
// ============================================================================
console.log("\n=== Application 4: Order Finding (Shor's Algorithm Preview) ===\n");

console.log("Shor's quantum algorithm for factoring n relies on finding the order");
console.log("of a random element a mod n. If ord_n(a) = r is even and a^(r/2) != -1,");
console.log("then gcd(a^(r/2) +/- 1, n) often gives a factor of n.\n");

function shorsMethod(n: bigint, a: bigint): { factor: bigint | null; steps: string[] } {
  const steps: string[] = [];

  if (gcd(a, n) > 1n) {
    steps.push(`Lucky! gcd(${a}, ${n}) = ${gcd(a, n)} is already a factor`);
    return { factor: gcd(a, n), steps };
  }

  const r = multiplicativeOrder(a, n);
  if (r === null) {
    steps.push(`${a} is not coprime to ${n}`);
    return { factor: null, steps };
  }

  steps.push(`Order of ${a} mod ${n}: r = ${r}`);

  if (r % 2n !== 0n) {
    steps.push(`r = ${r} is odd, method fails`);
    return { factor: null, steps };
  }

  const halfPower = power_mod(a, r / 2n, n);
  steps.push(`a^(r/2) mod n = ${a}^${r / 2n} mod ${n} = ${halfPower}`);

  if (halfPower === n - 1n) {
    steps.push(`a^(r/2) = -1 mod n, method fails`);
    return { factor: null, steps };
  }

  const factor1 = gcd(halfPower + 1n, n);
  const factor2 = gcd(halfPower - 1n, n);
  steps.push(`gcd(${halfPower} + 1, ${n}) = ${factor1}`);
  steps.push(`gcd(${halfPower} - 1, ${n}) = ${factor2}`);

  const factor = (factor1 > 1n && factor1 < n) ? factor1 :
                 (factor2 > 1n && factor2 < n) ? factor2 : null;

  return { factor, steps };
}

const compositeN = 15n; // 3 * 5
console.log(`Factoring ${compositeN} = ${formatFactorization(factor(compositeN))}:\n`);

for (let a = 2n; a < compositeN; a++) {
  if (gcd(a, compositeN) === 1n) {
    const result = shorsMethod(compositeN, a);
    console.log(`a = ${a}:`);
    for (const step of result.steps) {
      console.log(`  ${step}`);
    }
    if (result.factor) {
      console.log(`  Found factor: ${result.factor}`);
    }
    console.log();
  }
}

// ============================================================================
// Application 5: Cyclic Structure of (Z/nZ)*
// ============================================================================
console.log("\n=== Application 5: Structure of (Z/nZ)* ===\n");

console.log("Knowing phi(n) helps understand the group structure of (Z/nZ)*.\n");

function analyzeGroup(n: bigint): void {
  const phi = euler_phi(n);
  const Zn = Zmod(n);
  const units = Zn.units();

  console.log(`(Z/${n}Z)*:`);
  console.log(`  Order: phi(${n}) = ${phi}`);
  console.log(`  Elements: {${units.map(u => u.value).join(", ")}}`);

  // Find all element orders
  const orderCounts: Map<bigint, number> = new Map();
  let primitiveRoots: bigint[] = [];

  for (const unit of units) {
    const order = multiplicativeOrder(unit.value, n)!;
    orderCounts.set(order, (orderCounts.get(order) || 0) + 1);
    if (order === phi) {
      primitiveRoots.push(unit.value);
    }
  }

  console.log(`  Orders present: {${Array.from(orderCounts.keys()).sort((a, b) => a < b ? -1 : 1).join(", ")}}`);

  if (primitiveRoots.length > 0) {
    console.log(`  Primitive roots: {${primitiveRoots.join(", ")}}`);
    console.log(`  The group is CYCLIC (has primitive roots)`);
  } else {
    console.log(`  NO primitive roots exist`);
    console.log(`  The group is NOT cyclic`);
  }
  console.log();
}

// Analyze various groups
analyzeGroup(7n);   // Cyclic (prime)
analyzeGroup(8n);   // Not cyclic
analyzeGroup(9n);   // Cyclic (prime power, p odd)
analyzeGroup(15n);  // Not cyclic (product of odd primes)

console.log("Theorem: (Z/nZ)* is cyclic iff n = 1, 2, 4, p^k, or 2p^k (p odd prime).");

// ============================================================================
// Application 6: Carmichael Function in RSA
// ============================================================================
console.log("\n=== Application 6: Carmichael Function in RSA ===\n");

console.log("For RSA, the Carmichael function lambda(n) can replace phi(n):");
console.log("  lambda(pq) = lcm(p-1, q-1)  (often smaller than phi(n))");
console.log("This gives a smaller private exponent d.\n");

function carmichael(nn: bigint): bigint {
  const factors = factor(nn);
  let result = 1n;

  for (const [p, e] of factors) {
    if (p === -1n) continue;

    let lambdaPe: bigint;
    if (p === 2n && e >= 3n) {
      lambdaPe = p ** (e - 2n);
    } else {
      lambdaPe = p ** (e - 1n) * (p - 1n);
    }

    result = (result * lambdaPe) / gcd(result, lambdaPe);
  }

  return result;
}

// Compare RSA with phi vs lambda
const pRsa = 61n;
const qRsa = 53n;
const nRsa = pRsa * qRsa;
const phiRsa = euler_phi(nRsa);
const lambdaRsa = carmichael(nRsa);
const eRsa = 17n;
const dPhi = inverse_mod(eRsa, phiRsa);
const dLambda = inverse_mod(eRsa, lambdaRsa);

console.log(`RSA with n = ${pRsa} * ${qRsa} = ${nRsa}:`);
console.log(`  phi(n) = ${phiRsa}`);
console.log(`  lambda(n) = ${lambdaRsa}`);
console.log(`  Ratio: phi/lambda = ${phiRsa / lambdaRsa}`);
console.log(`\n  With e = ${eRsa}:`);
console.log(`    d using phi:    ${dPhi}`);
console.log(`    d using lambda: ${dLambda}`);

// Verify both work
const msgRsa = 42n;
const ctRsa = power_mod(msgRsa, eRsa, nRsa);
const decPhi = power_mod(ctRsa, dPhi, nRsa);
const decLambda = power_mod(ctRsa, dLambda, nRsa);

console.log(`\n  Decryption test (m = ${msgRsa}):`);
console.log(`    Using d from phi:    ${decPhi} [${decPhi === msgRsa ? "OK" : "FAIL"}]`);
console.log(`    Using d from lambda: ${decLambda} [${decLambda === msgRsa ? "OK" : "FAIL"}]`);

console.log(`\nUsing lambda gives a smaller private exponent,`);
console.log(`which can speed up decryption.`);

// ============================================================================
// Summary
// ============================================================================
console.log("\n" + "=".repeat(70));
console.log("SUMMARY");
console.log("=".repeat(70));
console.log(`
Applications of Euler's Theorem:

1. RSA Cryptosystem
   - Key generation uses phi(n) = (p-1)(q-1)
   - Correctness relies on a^phi(n) = 1 mod n
   - Can use lambda(n) for smaller private key

2. Multiplicative Order
   - ord_n(a) divides phi(n) (Lagrange's theorem)
   - Computed by testing divisors of phi(n)

3. Primitive Roots
   - Generator of (Z/nZ)* has order phi(n)
   - Number of primitive roots = phi(phi(n))
   - Exist iff n = 1, 2, 4, p^k, or 2p^k

4. Shor's Algorithm (Quantum)
   - Order finding is the key quantum subroutine
   - Classical part uses gcd(a^(r/2) +/- 1, n)

5. Group Structure Analysis
   - phi(n) gives the group order
   - Helps determine if group is cyclic

6. Carmichael Function
   - lambda(n) divides phi(n)
   - Universal exponent: a^lambda(n) = 1 for all units
   - Optimization for RSA private key
`);
