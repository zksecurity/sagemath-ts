/**
 * # Euler's Theorem
 *
 * Euler's theorem is a generalization of Fermat's Little Theorem to
 * arbitrary moduli (not just primes).
 *
 * ## Statement
 *
 * If gcd(a, n) = 1, then:
 *
 *   a^phi(n) = 1 (mod n)
 *
 * where phi(n) is Euler's totient function.
 *
 * ## Special Case: Fermat's Little Theorem
 *
 * When n = p is prime, phi(p) = p - 1, so:
 *   a^(p-1) = 1 (mod p)
 *
 * ## Proof Intuition
 *
 * The multiplicative group (Z/nZ)* has order phi(n).
 * By Lagrange's theorem, every element's order divides the group order.
 * Therefore a^phi(n) = 1 for all a in (Z/nZ)*.
 *
 * @module tutorial/part2/eulers-theorem
 */

import {
  euler_phi,
  power_mod,
  gcd,
  factor,
  Zmod,
  formatFactorization,
} from 'sagemath-ts';

console.log("=".repeat(70));
console.log("EULER'S THEOREM");
console.log("=".repeat(70));

// ============================================================================
// Example 1: Verifying Euler's theorem for various moduli
// ============================================================================
console.log("\n=== Example 1: Verifying a^phi(n) = 1 (mod n) ===\n");

const moduli = [8n, 12n, 15n, 20n];

for (const n of moduli) {
  const phi = euler_phi(n);
  console.log(`n = ${n}, phi(${n}) = ${phi}`);
  console.log(`Testing a^${phi} mod ${n} for all a coprime to ${n}:`);

  const results: string[] = [];
  for (let a = 1n; a < n; a++) {
    if (gcd(a, n) === 1n) {
      const result = power_mod(a, phi, n);
      results.push(`${a}^${phi}=${result}`);
    }
  }
  console.log(`  ${results.join(", ")}`);
  console.log(`  All equal to 1: ${results.every(r => r.endsWith("=1"))}\n`);
}

// ============================================================================
// Example 2: Comparison with Fermat's Little Theorem
// ============================================================================
console.log("\n=== Example 2: Euler's Theorem Generalizes Fermat ===\n");

console.log("For prime p, phi(p) = p - 1, so Euler's theorem becomes Fermat's theorem.\n");

const primes = [5n, 7n, 11n];

for (const p of primes) {
  const phi = euler_phi(p);
  console.log(`p = ${p}:`);
  console.log(`  phi(${p}) = ${phi} = ${p} - 1 (confirming p is prime)`);
  console.log(`  Euler: a^${phi} = 1 (mod ${p})`);
  console.log(`  Fermat: a^${p - 1n} = 1 (mod ${p})`);
  console.log(`  Same statement!`);

  // Verify with a = 2
  const a = 2n;
  console.log(`  Example: ${a}^${phi} mod ${p} = ${power_mod(a, phi, p)}\n`);
}

// ============================================================================
// Example 3: Why gcd(a, n) = 1 is required
// ============================================================================
console.log("\n=== Example 3: Why gcd(a, n) = 1 is Required ===\n");

const n3 = 12n;
const phi3 = euler_phi(n3);

console.log(`n = ${n3}, phi(${n3}) = ${phi3}`);
console.log(`\nTesting all a from 1 to ${n3 - 1n}:`);
console.log("-".repeat(60));

for (let a = 1n; a < n3; a++) {
  const g = gcd(a, n3);
  const result = power_mod(a, phi3, n3);
  const coprime = g === 1n;
  const equalsOne = result === 1n;
  const status = coprime ? (equalsOne ? "OK (coprime, = 1)" : "UNEXPECTED") :
                          (equalsOne ? "Coincidence" : "Expected != 1");
  console.log(`  a = ${a.toString().padStart(2)}: gcd = ${g}, ` +
              `${a}^${phi3} mod ${n3} = ${result.toString().padStart(2)}  [${status}]`);
}

console.log(`
When gcd(a, n) > 1, a is NOT in the multiplicative group (Z/nZ)*
and Euler's theorem does not apply.
`);

// ============================================================================
// Example 4: Order of elements and Euler's theorem
// ============================================================================
console.log("\n=== Example 4: Element Orders Divide phi(n) ===\n");

console.log("By Lagrange's theorem, the order of any element divides |G| = phi(n).");
console.log("Euler's theorem is a consequence: a^phi(n) = a^(order * k) = 1\n");

const n4 = 15n;
const phi4 = euler_phi(n4);
const Z15 = Zmod(n4);

console.log(`In (Z/${n4}Z)*, phi(${n4}) = ${phi4} = ${formatFactorization(factor(phi4))}`);
console.log(`Possible element orders must divide ${phi4}: {${[1n, 2n, 4n].join(", ")}}\n`);

console.log("Element | Order | Divides phi(n)?");
console.log("--------+-------+----------------");

for (const unit of Z15.units()) {
  const a = unit.value;
  let order = 1n;
  let current = a;
  while (current !== 1n) {
    current = (current * a) % n4;
    order++;
  }
  const divides = phi4 % order === 0n;
  console.log(`   ${a.toString().padStart(2)}    |   ${order}   | ${divides}`);
}

// ============================================================================
// Example 5: Computing modular inverses
// ============================================================================
console.log("\n=== Example 5: Modular Inverses via Euler's Theorem ===\n");

console.log("Since a^phi(n) = 1 (mod n), we have:");
console.log("  a * a^(phi(n) - 1) = 1 (mod n)");
console.log("  Therefore: a^(-1) = a^(phi(n) - 1) (mod n)\n");

const n5 = 20n;
const phi5 = euler_phi(n5);

console.log(`In Z/${n5}Z, phi(${n5}) = ${phi5}:`);
console.log(`Inverse of a is a^${phi5 - 1n} = a^${phi5 - 1n} mod ${n5}\n`);

console.log("  a | a^(-1) = a^" + (phi5 - 1n) + " | Verification");
console.log("----+--------+-------------");

for (let a = 1n; a < n5; a++) {
  if (gcd(a, n5) === 1n) {
    const inv = power_mod(a, phi5 - 1n, n5);
    const verify = (a * inv) % n5;
    console.log(` ${a.toString().padStart(2)} |   ${inv.toString().padStart(2)}   | ${a} * ${inv} mod ${n5} = ${verify}`);
  }
}

// ============================================================================
// Example 6: Exponent reduction with Euler's theorem
// ============================================================================
console.log("\n=== Example 6: Exponent Reduction ===\n");

console.log("When computing a^k mod n with gcd(a, n) = 1:");
console.log("  a^k = a^(k mod phi(n)) (mod n)");
console.log("This dramatically reduces computation for large exponents.\n");

const n6 = 21n; // 3 * 7
const phi6 = euler_phi(n6); // (3-1)(7-1) = 12
const a6 = 5n;
const bigExp = 1000000000000000000000n; // 10^21

console.log(`Computing ${a6}^${bigExp} mod ${n6}:`);
console.log(`  phi(${n6}) = ${phi6}`);
console.log(`  ${bigExp} mod ${phi6} = ${bigExp % phi6}`);
console.log(`  Therefore ${a6}^${bigExp} mod ${n6} = ${a6}^${bigExp % phi6} mod ${n6}`);
console.log(`  = ${power_mod(a6, bigExp % phi6, n6)}`);

// Verify (power_mod handles this internally)
const directResult = power_mod(a6, bigExp, n6);
console.log(`\n  Direct computation: ${directResult}`);
console.log(`  Results match: ${directResult === power_mod(a6, bigExp % phi6, n6)}`);

// ============================================================================
// Example 7: The Carmichael function (reduced exponent)
// ============================================================================
console.log("\n=== Example 7: Carmichael Function lambda(n) ===\n");

console.log("The Carmichael function lambda(n) is the smallest positive m such that:");
console.log("  a^m = 1 (mod n) for ALL a coprime to n");
console.log("lambda(n) divides phi(n), and can be strictly smaller.\n");

function carmichael(n: bigint): bigint {
  if (n === 1n) return 1n;
  if (n === 2n) return 1n;
  if (n === 4n) return 2n;

  const factors = factor(n);
  let result = 1n;

  for (const [p, e] of factors) {
    if (p === -1n) continue;

    let lambda_pe: bigint;
    if (p === 2n && e >= 3n) {
      // Special case for powers of 2 >= 8
      lambda_pe = p ** (e - 2n);
    } else {
      lambda_pe = p ** (e - 1n) * (p - 1n);
    }

    // LCM of result and lambda_pe
    result = (result * lambda_pe) / gcd(result, lambda_pe);
  }

  return result;
}

const carmichaelExamples = [8n, 12n, 15n, 20n, 24n, 30n];

console.log("   n | phi(n) | lambda(n) | ratio");
console.log("-----+--------+-----------+------");
for (const n of carmichaelExamples) {
  const phi = euler_phi(n);
  const lambda = carmichael(n);
  const ratio = phi / lambda;
  console.log(` ${n.toString().padStart(3)} | ${phi.toString().padStart(6)} | ${lambda.toString().padStart(9)} | ${ratio}`);
}

console.log("\nExample: Verifying lambda(24) = 2");
const n7 = 24n;
const lambda7 = carmichael(n7);
const Z24 = Zmod(n7);

console.log(`All units in Z/${n7}Z raised to power ${lambda7}:`);
for (const unit of Z24.units()) {
  const result = power_mod(unit.value, lambda7, n7);
  console.log(`  ${unit.value}^${lambda7} mod ${n7} = ${result}`);
}

// ============================================================================
// Summary
// ============================================================================
console.log("\n" + "=".repeat(70));
console.log("SUMMARY");
console.log("=".repeat(70));
console.log(`
Euler's Theorem:
  If gcd(a, n) = 1, then a^phi(n) = 1 (mod n)

Key Points:
1. Generalizes Fermat's Little Theorem to composite moduli
2. Requires a to be coprime to n (a must be a unit in Z/nZ)
3. phi(n) = |(Z/nZ)*| is the order of the multiplicative group
4. Element orders divide phi(n) (Lagrange's theorem)

Applications:
1. Modular inverse: a^(-1) = a^(phi(n)-1) mod n
2. Exponent reduction: a^k = a^(k mod phi(n)) mod n
3. RSA encryption: works because m^(ed) = m^(1 + k*phi(n)) = m

Carmichael Function:
- lambda(n) <= phi(n), can be strictly smaller
- Gives the actual universal exponent for (Z/nZ)*
- Important for RSA security analysis
`);
