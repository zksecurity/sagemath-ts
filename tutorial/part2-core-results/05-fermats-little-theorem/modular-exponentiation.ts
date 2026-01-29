/**
 * # Fast Modular Exponentiation
 *
 * Computing a^n mod m efficiently is fundamental to cryptography.
 * Naive exponentiation is O(n) multiplications, but we can do much better.
 *
 * ## Square-and-Multiply Algorithm
 *
 * The key insight: express n in binary and use the property:
 *   a^(2k) = (a^k)^2
 *   a^(2k+1) = a * (a^k)^2
 *
 * This reduces the computation to O(log n) multiplications.
 *
 * ## Connection to Fermat's Theorem
 *
 * When computing a^n mod p (for prime p), we can first reduce:
 *   n' = n mod (p-1)
 * Because a^(p-1) = 1 (mod p), so a^n = a^(n mod (p-1)) (mod p)
 *
 * This is especially useful when n is very large.
 *
 * @module tutorial/part2/modular-exponentiation
 */

import {
  power_mod,
  euler_phi,
  is_prime,
  Zmod,
} from 'sagemath-ts';

console.log("=".repeat(70));
console.log("FAST MODULAR EXPONENTIATION");
console.log("=".repeat(70));

// ============================================================================
// Example 1: The naive approach (for educational purposes)
// ============================================================================
console.log("\n=== Example 1: Naive vs. Fast Exponentiation ===\n");

function naiveModPow(a: bigint, n: bigint, m: bigint): bigint {
  let result = 1n;
  a = ((a % m) + m) % m;
  for (let i = 0n; i < n; i++) {
    result = (result * a) % m;
  }
  return result;
}

function squareAndMultiply(a: bigint, n: bigint, m: bigint): { result: bigint; steps: string[] } {
  const steps: string[] = [];
  let result = 1n;
  a = ((a % m) + m) % m;

  steps.push(`Starting: result = 1, base = ${a}`);
  steps.push(`Binary representation of ${n}: ${n.toString(2)}`);

  let bitPosition = 0;
  let base = a;
  while (n > 0n) {
    if ((n & 1n) === 1n) {
      result = (result * base) % m;
      steps.push(`Bit ${bitPosition} is 1: result = result * ${base} mod ${m} = ${result}`);
    } else {
      steps.push(`Bit ${bitPosition} is 0: result unchanged = ${result}`);
    }
    base = (base * base) % m;
    if (n > 1n) {
      steps.push(`             square base: ${base}`);
    }
    n >>= 1n;
    bitPosition++;
  }

  return { result, steps };
}

const a1 = 3n;
const n1 = 13n;
const m1 = 17n;

console.log(`Computing ${a1}^${n1} mod ${m1}:`);
console.log(`\nNaive approach: ${n1} multiplications`);

const { result, steps } = squareAndMultiply(a1, n1, m1);
console.log(`\nSquare-and-multiply approach:`);
for (const step of steps) {
  console.log(`  ${step}`);
}
console.log(`\nFinal result: ${result}`);
console.log(`Number of multiplications: ${steps.filter(s => s.includes("result =") || s.includes("square")).length}`);
console.log(`\nVerification with power_mod: ${power_mod(a1, n1, m1)}`);

// ============================================================================
// Example 2: Using Fermat's theorem to reduce the exponent
// ============================================================================
console.log("\n=== Example 2: Exponent Reduction with Fermat's Theorem ===\n");

const p = 17n; // prime
const base = 5n;
const hugeExponent = 10000000000000000000n;

console.log(`Computing ${base}^${hugeExponent} mod ${p}`);
console.log(`\nDirect computation would need ~${hugeExponent.toString().length} digit operations!`);

// Using Fermat's theorem: a^(p-1) = 1 (mod p)
// So a^n = a^(n mod (p-1)) (mod p)
const reducedExp = hugeExponent % (p - 1n);
console.log(`\nUsing Fermat's theorem:`);
console.log(`  Since p = ${p} is prime, ${base}^${p - 1n} = 1 (mod ${p})`);
console.log(`  Therefore ${base}^${hugeExponent} = ${base}^(${hugeExponent} mod ${p - 1n})`);
console.log(`                    = ${base}^${reducedExp} (mod ${p})`);

const result2 = power_mod(base, reducedExp, p);
console.log(`\n  Result: ${result2}`);

// Verify with direct computation using power_mod (which handles this efficiently)
const directResult = power_mod(base, hugeExponent, p);
console.log(`  Verification: ${directResult}`);
console.log(`  Match: ${result2 === directResult}`);

// ============================================================================
// Example 3: Computing modular inverses using Fermat's theorem
// ============================================================================
console.log("\n=== Example 3: Modular Inverses via Fermat's Theorem ===\n");

console.log("For prime p and a coprime to p:");
console.log("  a^(p-1) = 1 (mod p)");
console.log("  a * a^(p-2) = 1 (mod p)");
console.log("  Therefore: a^(-1) = a^(p-2) (mod p)\n");

const p3 = 13n;
console.log(`In Z/${p3}Z, computing inverses:`);
console.log("-".repeat(40));

for (let a = 1n; a < p3; a++) {
  const inv = power_mod(a, p3 - 2n, p3);
  const check = (a * inv) % p3;
  console.log(`  ${a}^(-1) = ${a}^${p3 - 2n} mod ${p3} = ${inv}  (verify: ${a} * ${inv} = ${check} mod ${p3})`);
}

// ============================================================================
// Example 4: Timing comparison
// ============================================================================
console.log("\n=== Example 4: Efficiency Demonstration ===\n");

// Compare different exponent sizes
const exponents = [
  { name: "Small", exp: 100n },
  { name: "Medium", exp: 1000000n },
  { name: "Large", exp: 1000000000000n },
  { name: "Huge", exp: 10n ** 50n },
  { name: "Massive", exp: 10n ** 100n },
];

const p4 = 1000000007n; // A large prime
const base4 = 12345n;

console.log(`Computing ${base4}^exp mod ${p4} for various exponents:\n`);
console.log("Exponent Size".padEnd(12) + "Exponent".padEnd(35) + "Result");
console.log("-".repeat(60));

for (const { name, exp } of exponents) {
  const expStr = exp.toString();
  const expDisplay = expStr.length > 30 ? `10^${(expStr.length - 1)}` : expStr;
  const result4 = power_mod(base4, exp, p4);
  console.log(`${name.padEnd(12)}${expDisplay.padEnd(35)}${result4}`);
}

console.log("\nNote: All computations complete instantly thanks to:");
console.log("  1. Square-and-multiply algorithm: O(log exp) multiplications");
console.log("  2. Exponent reduction: exp mod phi(p) before computation");

// ============================================================================
// Example 5: Working with the ring structure
// ============================================================================
console.log("\n=== Example 5: Exponentiation with Ring Elements ===\n");

const Z23 = Zmod(23n);
const elem = Z23.__call__(7n);

console.log(`Working in Z/23Z with element ${elem.value}:`);
console.log(`\nPowers of ${elem.value}:`);

let current = Z23.__call__(1n);
const powers: bigint[] = [];
for (let i = 0; i <= 22; i++) {
  powers.push(current.value);
  current = current.mul(elem);
}

// Display in rows of 8
for (let i = 0; i <= 22; i += 8) {
  const row = powers.slice(i, Math.min(i + 8, 23))
    .map((v, j) => `7^${(i + j).toString().padStart(2)} = ${v.toString().padStart(2)}`)
    .join("  ");
  console.log(`  ${row}`);
}

console.log(`\nNote: 7^22 = ${powers[22]} = 1 (Fermat's theorem)`);
console.log(`The order of 7 in (Z/23Z)* is 22 (7 is a primitive root)`);

// ============================================================================
// Example 6: Cryptographic applications preview
// ============================================================================
console.log("\n=== Example 6: Cryptographic Application Preview ===\n");

console.log("RSA Key Generation (simplified example):");
console.log("-".repeat(50));

// Small primes for demonstration
const pp = 61n;
const qq = 53n;
const nn = pp * qq;
const phiN = (pp - 1n) * (qq - 1n);
const ee = 17n; // Public exponent
// Private exponent d such that e*d = 1 (mod phi(n))
const dd = power_mod(ee, phiN - 1n, phiN); // Using Fermat since phi(n) might not be prime

console.log(`p = ${pp}, q = ${qq}`);
console.log(`n = p * q = ${nn}`);
console.log(`phi(n) = (p-1)(q-1) = ${phiN}`);
console.log(`e = ${ee} (public exponent)`);
console.log(`d = e^(-1) mod phi(n) = ${dd} (private exponent)`);

// Encryption and decryption
const message = 42n;
const ciphertext = power_mod(message, ee, nn);
const decrypted = power_mod(ciphertext, dd, nn);

console.log(`\nEncryption/Decryption:`);
console.log(`  Message m = ${message}`);
console.log(`  Ciphertext c = m^e mod n = ${message}^${ee} mod ${nn} = ${ciphertext}`);
console.log(`  Decrypted m' = c^d mod n = ${ciphertext}^${dd} mod ${nn} = ${decrypted}`);
console.log(`  Correct: ${message === decrypted}`);

console.log("\nWhy this works (using Fermat/Euler):");
console.log(`  c^d = (m^e)^d = m^(ed) mod n`);
console.log(`  Since ed = 1 (mod phi(n)), ed = 1 + k*phi(n) for some k`);
console.log(`  m^(ed) = m^(1 + k*phi(n)) = m * (m^phi(n))^k`);
console.log(`  By Euler's theorem, m^phi(n) = 1 (mod n)`);
console.log(`  Therefore m^(ed) = m * 1^k = m (mod n)`);

// ============================================================================
// Summary
// ============================================================================
console.log("\n" + "=".repeat(70));
console.log("SUMMARY");
console.log("=".repeat(70));
console.log(`
Key Points:
1. Square-and-multiply reduces O(n) to O(log n) multiplications
2. Fermat's theorem allows exponent reduction: a^n = a^(n mod (p-1)) (mod p)
3. Modular inverse via Fermat: a^(-1) = a^(p-2) (mod p)
4. These techniques are essential for cryptographic operations
5. The power_mod function in sagemath-ts handles all optimizations

Algorithm Complexity:
- Naive: O(n) multiplications
- Square-and-multiply: O(log n) multiplications
- Each multiplication of k-bit numbers: O(k^2) or O(k log k) with FFT
- Total for k-bit modulus and n-bit exponent: O(n * k^2) or O(n * k log k)
`);
