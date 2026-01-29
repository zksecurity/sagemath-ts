/**
 * # Modular Operations: Addition, Subtraction, and Multiplication
 *
 * Modular arithmetic operations form the computational backbone of cryptography.
 * Every RSA encryption, every elliptic curve point addition, and every
 * hash function computation relies on efficient modular operations.
 *
 * ## Mathematical Background
 *
 * For integers a, b and positive modulus n:
 *   - (a + b) mod n = ((a mod n) + (b mod n)) mod n
 *   - (a - b) mod n = ((a mod n) - (b mod n)) mod n
 *   - (a * b) mod n = ((a mod n) * (b mod n)) mod n
 *   - a^k mod n can be computed efficiently using square-and-multiply
 *
 * The key insight: we can reduce intermediate results modulo n at any point
 * without affecting the final answer. This keeps numbers small!
 *
 * ## Why This Matters for Cryptography
 *
 *   - RSA: c = m^e mod n (encryption), m = c^d mod n (decryption)
 *   - Diffie-Hellman: g^x mod p
 *   - ECDSA: Computations in a prime field
 *   - AES: Operations in GF(2^8)
 */

import { gcd, power_mod, Mod, Zmod } from 'sagemath-ts';

// =============================================================================
// Example 1: Basic Modular Addition and Subtraction
// =============================================================================
console.log("=== Example 1: Modular Addition and Subtraction ===\n");

/**
 * Modular addition and subtraction work like a clock.
 * If you go past n-1, you wrap around to 0.
 * If you go below 0, you wrap around to n-1.
 */
const Z12 = Zmod(12);  // Like a 12-hour clock

console.log("Clock arithmetic in Z/12Z:\n");

// Addition
const hour1 = Z12.__call__(9);
const hour2 = Z12.__call__(5);
const sumHours = hour1.add(hour2);
console.log(`  9 + 5 = ${sumHours.value} (mod 12)  [9 o'clock + 5 hours = ${sumHours.value} o'clock]`);

// Subtraction
const diff = hour1.sub(hour2);
console.log(`  9 - 5 = ${diff.value} (mod 12)`);

// Subtraction with wrap-around
const hour3 = Z12.__call__(3);
const diff2 = hour3.sub(hour2);
console.log(`  3 - 5 = ${diff2.value} (mod 12)  [going 5 hours before 3 o'clock = ${diff2.value} o'clock]`);

// Chain of operations
const result = Z12.__call__(7).add(8).sub(3);
console.log(`  7 + 8 - 3 = ${result.value} (mod 12)`);

console.log();

// =============================================================================
// Example 2: Modular Multiplication
// =============================================================================
console.log("=== Example 2: Modular Multiplication ===\n");

/**
 * Modular multiplication is the workhorse of cryptography.
 * The key property: (a * b) mod n = ((a mod n) * (b mod n)) mod n
 */
const Z17 = Zmod(17);

console.log("Multiplication in Z/17Z:\n");

const a = Z17.__call__(7);
const b = Z17.__call__(5);
const product = a.mul(b);

console.log(`  ${a.value} * ${b.value} = ${7n * 5n} = ${product.value} (mod 17)`);

// Building a multiplication table
console.log("\nMultiplication table for Z/5Z:");
const Z5 = Zmod(5);
console.log("    |  0  1  2  3  4");
console.log("----+---------------");
for (let i = 0n; i < 5n; i++) {
  let row = `  ${i} |`;
  for (let j = 0n; j < 5n; j++) {
    const prod = Z5.__call__(i).mul(j);
    row += `  ${prod.value}`;
  }
  console.log(row);
}

console.log();

// =============================================================================
// Example 3: Keeping Numbers Small
// =============================================================================
console.log("=== Example 3: Reducing Intermediate Results ===\n");

/**
 * A crucial optimization: reduce after each operation to prevent overflow.
 * This is how cryptographic libraries handle large numbers efficiently.
 */
const modulus = 1000000007n;  // A common prime in competitive programming

console.log(`Computing products modulo ${modulus}:\n`);

// Bad way: compute full product first, then reduce
// (This would overflow for large numbers!)
function badMethod(a: bigint, b: bigint, n: bigint): bigint {
  return (a * b) % n;  // Fine for small numbers, but can overflow!
}

// Good way: reduce inputs first
function goodMethod(a: bigint, b: bigint, n: bigint): bigint {
  return ((a % n) * (b % n)) % n;
}

const big1 = 123456789012345n;
const big2 = 987654321098765n;

console.log(`  a = ${big1}`);
console.log(`  b = ${big2}`);
console.log(`  a mod n = ${big1 % modulus}`);
console.log(`  b mod n = ${big2 % modulus}`);
console.log(`  (a * b) mod n = ${goodMethod(big1, big2, modulus)}`);

console.log();

// =============================================================================
// Example 4: Modular Exponentiation (Square-and-Multiply)
// =============================================================================
console.log("=== Example 4: Fast Modular Exponentiation ===\n");

/**
 * Computing a^n mod m directly would be impossibly slow for cryptographic sizes.
 * The square-and-multiply algorithm computes this in O(log n) multiplications.
 *
 * Key idea: a^13 = a^(1101 binary) = a^8 * a^4 * a^1
 */

// Trace through the square-and-multiply algorithm
function powerModTrace(base: bigint, exp: bigint, mod: bigint): bigint {
  console.log(`Computing ${base}^${exp} mod ${mod}:`);
  console.log(`  Binary exponent: ${exp.toString(2)}`);

  let result = 1n;
  let currentPower = base % mod;
  let tempExp = exp;
  let bitPosition = 0;

  console.log("\n  Step-by-step:");
  while (tempExp > 0n) {
    if ((tempExp & 1n) === 1n) {
      const oldResult = result;
      result = (result * currentPower) % mod;
      console.log(`    Bit ${bitPosition} is 1: result = ${oldResult} * ${currentPower} = ${result} (mod ${mod})`);
    } else {
      console.log(`    Bit ${bitPosition} is 0: skip`);
    }

    const oldPower = currentPower;
    currentPower = (currentPower * currentPower) % mod;
    console.log(`    Square: ${oldPower}^2 = ${currentPower} (mod ${mod})`);

    tempExp >>= 1n;
    bitPosition++;
  }

  console.log(`\n  Result: ${result}`);
  return result;
}

powerModTrace(3n, 13n, 17n);

console.log("\nUsing library function:");
console.log(`  power_mod(3, 13, 17) = ${power_mod(3n, 13n, 17n)}`);

console.log();

// =============================================================================
// Example 5: Large Exponentiation
// =============================================================================
console.log("=== Example 5: Cryptographic-Scale Exponentiation ===\n");

/**
 * In real cryptography, exponents can be hundreds of bits long.
 * Square-and-multiply handles this efficiently.
 */
const cryptoMod = 2n ** 127n - 1n;  // A Mersenne prime
const base = 12345n;
const largeExp = 2n ** 100n - 1n;  // A 100-bit number!

console.log(`Computing ${base}^(2^100 - 1) mod (2^127 - 1):`);
console.log(`  Exponent has ${largeExp.toString(2).length} bits`);

const start = performance.now();
const powerResult = power_mod(base, largeExp, cryptoMod);
const elapsed = performance.now() - start;

console.log(`  Result: ${powerResult}`);
console.log(`  Computed in: ${elapsed.toFixed(2)} ms`);
console.log("\nWithout square-and-multiply, this would require 2^100 multiplications!");

console.log();

// =============================================================================
// Example 6: Modular Arithmetic with IntegerMod
// =============================================================================
console.log("=== Example 6: Using the IntegerMod Class ===\n");

/**
 * The IntegerMod class provides a convenient interface for modular arithmetic.
 * All operations automatically reduce results modulo n.
 */
const p = 23n;
const Zp = Zmod(p);

console.log(`Working in Z/${p}Z:\n`);

const x = Zp.__call__(7);
const y = Zp.__call__(15);

console.log(`  x = ${x.value}`);
console.log(`  y = ${y.value}`);
console.log(`  x + y = ${x.add(y).value}`);
console.log(`  x - y = ${x.sub(y).value}`);
console.log(`  x * y = ${x.mul(y).value}`);
console.log(`  x^5 = ${x.pow(5).value}`);

// Chaining operations
const chain = Zp.__call__(2).mul(3).add(5).pow(2);
console.log(`\n  (2 * 3 + 5)^2 = ${chain.value} (mod 23)`);

console.log();

// =============================================================================
// Example 7: Distributive Law in Modular Arithmetic
// =============================================================================
console.log("=== Example 7: Algebraic Properties ===\n");

/**
 * Modular arithmetic preserves the familiar algebraic properties:
 *   - Commutative: a + b = b + a, a * b = b * a
 *   - Associative: (a + b) + c = a + (b + c)
 *   - Distributive: a * (b + c) = a*b + a*c
 */
const Z13 = Zmod(13);

const vals = [3n, 7n, 11n];
const [va, vb, vc] = vals.map(v => Z13.__call__(v));

console.log(`Testing algebraic properties in Z/13Z with a=${vals[0]}, b=${vals[1]}, c=${vals[2]}:\n`);

// Commutative (addition)
console.log(`  Commutative (add): ${va.add(vb).value} = ${vb.add(va).value} ? ${va.add(vb).value === vb.add(va).value}`);

// Commutative (multiplication)
console.log(`  Commutative (mul): ${va.mul(vb).value} = ${vb.mul(va).value} ? ${va.mul(vb).value === vb.mul(va).value}`);

// Associative (addition)
const leftAssocAdd = va.add(vb).add(vc);
const rightAssocAdd = va.add(vb.add(vc));
console.log(`  Associative (add): (a+b)+c = ${leftAssocAdd.value}, a+(b+c) = ${rightAssocAdd.value} ? ${leftAssocAdd.value === rightAssocAdd.value}`);

// Distributive
const distLeft = va.mul(vb.add(vc));
const distRight = va.mul(vb).add(va.mul(vc));
console.log(`  Distributive: a*(b+c) = ${distLeft.value}, a*b + a*c = ${distRight.value} ? ${distLeft.value === distRight.value}`);

console.log();

// =============================================================================
// Example 8: Negative Exponents
// =============================================================================
console.log("=== Example 8: Negative Exponents ===\n");

/**
 * a^(-k) mod n = (a^(-1))^k mod n, where a^(-1) is the modular inverse.
 * This only works when gcd(a, n) = 1.
 */
const Z29 = Zmod(29);  // Prime modulus, so all non-zero elements are invertible

const elem = Z29.__call__(5);

console.log(`Working with 5 in Z/29Z:\n`);
console.log(`  5^1 = ${elem.pow(1).value}`);
console.log(`  5^2 = ${elem.pow(2).value}`);
console.log(`  5^(-1) = ${elem.pow(-1).value}`);
console.log(`  5^(-2) = ${elem.pow(-2).value}`);

// Verify: 5^2 * 5^(-2) = 1
const verify = elem.pow(2).mul(elem.pow(-2));
console.log(`\n  Verification: 5^2 * 5^(-2) = ${verify.value} (should be 1)`);

console.log();

// =============================================================================
// Example 9: Modular Arithmetic Performance
// =============================================================================
console.log("=== Example 9: Performance Considerations ===\n");

/**
 * For cryptographic applications, we need efficient implementations.
 * Key optimizations:
 *   1. Reduce after each operation
 *   2. Use square-and-multiply for exponentiation
 *   3. Use Montgomery multiplication for repeated operations
 */
console.log("Comparing different exponentiation strategies:\n");

const testMod = 2n ** 61n - 1n;  // Mersenne prime
const testBase = 123456789n;
const testExp = 10000000n;

// Method 1: Library function (optimized)
const start1 = performance.now();
const result1 = power_mod(testBase, testExp, testMod);
const time1 = performance.now() - start1;

console.log(`  power_mod(${testBase}, ${testExp}, ${testMod}):`);
console.log(`    Result: ${result1}`);
console.log(`    Time: ${time1.toFixed(2)} ms`);

// Method 2: Using Mod class
const start2 = performance.now();
const result2 = Mod(testBase, testMod).pow(testExp);
const time2 = performance.now() - start2;

console.log(`\n  Mod(${testBase}, ${testMod}).pow(${testExp}):`);
console.log(`    Result: ${result2.value}`);
console.log(`    Time: ${time2.toFixed(2)} ms`);

console.log();

// =============================================================================
// Exercise: Implement Modular Exponentiation
// =============================================================================
console.log("=== Exercise: Implement Square-and-Multiply ===\n");

/**
 * EXERCISE: Implement the square-and-multiply algorithm for modular exponentiation.
 * This is one of the most important algorithms in cryptography!
 */

// TODO: Implement this function
function myPowerMod(base: bigint, exp: bigint, mod: bigint): bigint {
  if (mod === 1n) return 0n;

  let result = 1n;
  base = base % mod;

  while (exp > 0n) {
    // If current bit is 1, multiply result by base
    if ((exp & 1n) === 1n) {
      result = (result * base) % mod;
    }
    // Square the base
    base = (base * base) % mod;
    // Move to next bit
    exp >>= 1n;
  }

  return result;
}

console.log("Testing your modular exponentiation:");
const testCases: [bigint, bigint, bigint][] = [
  [3n, 13n, 17n],
  [2n, 10n, 1000n],
  [7n, 100n, 13n],
];

for (const [b, e, m] of testCases) {
  const yours = myPowerMod(b, e, m);
  const library = power_mod(b, e, m);
  console.log(`  ${b}^${e} mod ${m}: yours=${yours}, library=${library}, match=${yours === library}`);
}

console.log();

// =============================================================================
// Cryptographic Application: RSA Encryption
// =============================================================================
console.log("=== Cryptographic Application: RSA Encryption ===\n");

/**
 * RSA encryption: ciphertext = message^e mod n
 * RSA decryption: message = ciphertext^d mod n
 *
 * All the modular arithmetic we've learned comes together here!
 */
// Small RSA parameters for demonstration
const rsaP = 61n;
const rsaQ = 53n;
const rsaN = rsaP * rsaQ;  // 3233
const phiN = (rsaP - 1n) * (rsaQ - 1n);  // 3120
const rsaE = 17n;  // Public exponent

// Compute private exponent d
let rsaD = 1n;
while ((rsaD * rsaE) % phiN !== 1n) rsaD++;

console.log("RSA Parameters:");
console.log(`  p = ${rsaP}, q = ${rsaQ}`);
console.log(`  n = ${rsaN}`);
console.log(`  e = ${rsaE} (public)`);
console.log(`  d = ${rsaD} (private)`);

// Encrypt a message
const message = 42n;
const ciphertext = power_mod(message, rsaE, rsaN);
const decrypted = power_mod(ciphertext, rsaD, rsaN);

console.log(`\nEncryption/Decryption:`);
console.log(`  Message: ${message}`);
console.log(`  Encrypted: ${message}^${rsaE} mod ${rsaN} = ${ciphertext}`);
console.log(`  Decrypted: ${ciphertext}^${rsaD} mod ${rsaN} = ${decrypted}`);
console.log(`\nRSA works because of Euler's theorem: a^phi(n) = 1 (mod n)`);
