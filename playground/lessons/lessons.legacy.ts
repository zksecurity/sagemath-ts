export type NavItem = {
  label: string;
  href?: string;
  enabled: boolean;
};

export type NavGroup = {
  title: string;
  items: NavItem[];
};

export type LessonBlock =
  | { type: 'paragraph'; text: string }
  | { type: 'code'; code: string };

export type Lesson = {
  id: string;
  title: string;
  blocks: LessonBlock[];
};

export const nav: NavGroup[] = [
  {
    title: 'Part 1: Foundations',
    items: [
      { label: '01. Integers and divisibility', href: '#lesson-gcd', enabled: true },
      { label: '02. Modular arithmetic', href: '#lesson-modinv', enabled: true },
      { label: '03. Groups', enabled: false },
      { label: '04. Rings and fields', enabled: false },
    ],
  },
  {
    title: 'Part 2: Core results',
    items: [
      { label: "05. Fermat's little theorem", enabled: false },
      { label: "06. Euler's theorem", enabled: false },
      { label: '07. Chinese remainder theorem', href: '#lesson-crt', enabled: true },
      { label: '08. Quadratic residues', href: '#lesson-qr', enabled: true },
    ],
  },
  {
    title: 'Part 3: Computational',
    items: [
      { label: '09. Primality testing', href: '#lesson-primes', enabled: true },
      { label: '10. Factorization', enabled: false },
      { label: '11. Discrete logarithm', enabled: false },
    ],
  },
  {
    title: 'Part 4: Public key',
    items: [
      { label: '12. RSA', href: '#lesson-rsa', enabled: true },
      { label: '13. Diffie-Hellman', enabled: false },
    ],
  },
  {
    title: 'Part 5: Elliptic curves',
    items: [
      { label: '14. Elliptic curves intro', href: '#lesson-ec', enabled: true },
      { label: '15. EC over finite fields', enabled: false },
      { label: '16. ECDH and ECDSA', enabled: false },
    ],
  },
  {
    title: 'Part 6: Advanced',
    items: [
      { label: '17. Pairings', enabled: false },
      { label: '18. Lattices intro', enabled: false },
    ],
  },
];

export const lessons: Lesson[] = [
  {
    id: 'gcd',
    title: '1. GCD and the Euclidean algorithm',
    blocks: [
      {
        type: 'paragraph',
        text:
          'The greatest common divisor of a and b is the largest integer dividing both. The Euclidean algorithm computes it using repeated remainders.',
      },
      {
        type: 'paragraph',
        text:
          'It is the backbone of modular inverses, RSA key generation, and many factoring methods.',
      },
      {
        type: 'paragraph',
        text:
          'The extended Euclidean algorithm also produces coefficients x and y such that ax + by = gcd(a, b). Those coefficients are what give us modular inverses.',
      },
      {
        type: 'code',
        code: `// Compute gcd for a few pairs\nconst { gcd } = Sage;\n\nconst pairs = [\n  [48n, 18n],\n  [252n, 198n],\n  [1071n, 462n],\n];\n\nfor (const [a, b] of pairs) {\n  print('gcd(' + a + ', ' + b + ') =', gcd(a, b));\n}`,
      },
    ],
  },
  {
    id: 'modinv',
    title: '2. Modular inverses',
    blocks: [
      {
        type: 'paragraph',
        text:
          'A modular inverse of a modulo n is a number x such that a x is congruent to 1 (mod n). It exists only if gcd(a, n) = 1.',
      },
      {
        type: 'paragraph',
        text:
          'Computationally, we use the extended Euclidean algorithm to solve ax + ny = 1, then reduce x into the range [0, n).',
      },
      {
        type: 'code',
        code: `// Compute a modular inverse and verify it\nconst { inverse_mod } = Sage;\n\nconst a = 37n;\nconst n = 101n;\nconst inv = inverse_mod(a, n);\n\nprint('inverse =', inv);\nprint('check =', (a * inv) % n);`,
      },
    ],
  },
  {
    id: 'crt',
    title: '3. Chinese remainder theorem',
    blocks: [
      {
        type: 'paragraph',
        text:
          'CRT combines several congruences into a single solution modulo the product. It is the backbone of many speedups in crypto.',
      },
      {
        type: 'paragraph',
        text:
          'When the moduli are pairwise coprime, there is a unique solution modulo their product. RSA decryption often uses CRT to speed up exponentiation.',
      },
      {
        type: 'code',
        code: `// Solve a classic CRT system\nconst { CRT_list } = Sage;\n\nconst residues = [2n, 3n, 2n];\nconst moduli = [3n, 5n, 7n];\nconst x = CRT_list(residues, moduli);\n\nprint('x =', x);\nprint('check:', x % 3n, x % 5n, x % 7n);`,
      },
    ],
  },
  {
    id: 'qr',
    title: '4. Quadratic residues',
    blocks: [
      {
        type: 'paragraph',
        text:
          'Quadratic residues tell us whether a has a square root mod p. They drive quadratic reciprocity, primality proofs, and EC arithmetic.',
      },
      {
        type: 'paragraph',
        text:
          'The Legendre symbol returns 1 if a is a residue, -1 if it is a non-residue, and 0 if a is divisible by p. sqrt_mod uses Tonelli-Shanks under the hood.',
      },
      {
        type: 'code',
        code: `// Legendre symbol demo\nconst { legendre_symbol, sqrt_mod } = Sage;\n\nconst p = 97n;\nconst a = 56n;\nconst symbol = legendre_symbol(a, p);\n\nprint('legendre_symbol(a, p) =', symbol);\nif (symbol === 1n) {\n  const r = sqrt_mod(a, p);\n  print('sqrt =', r, 'check =', (r * r) % p);\n} else {\n  print('no square root for a mod p');\n}`,
      },
    ],
  },
  {
    id: 'primes',
    title: '5. Primality testing',
    blocks: [
      {
        type: 'paragraph',
        text:
          'Probabilistic tests like Miller-Rabin let us test huge numbers quickly with tiny error probability.',
      },
      {
        type: 'paragraph',
        text:
          'A small set of bases makes Miller-Rabin deterministic for large ranges. sagemath-ts uses the deterministic bases from Sage for safe results.',
      },
      {
        type: 'code',
        code: `// Compare a few numbers with is_prime\nconst { is_prime } = Sage;\n\nconst numbers = [97n, 561n, 1105n, 1009n, 1013n];\nfor (const n of numbers) {\n  print(n, 'is prime?', is_prime(n));\n}`,
      },
    ],
  },
  {
    id: 'rsa',
    title: '6. RSA round-trip',
    blocks: [
      {
        type: 'paragraph',
        text:
          'RSA uses modular exponentiation in the multiplicative group of a modulus n = pq. This example uses toy values.',
      },
      {
        type: 'paragraph',
        text:
          'Real RSA requires large primes and padding. This example is only for understanding the flow of key generation, encryption, and decryption.',
      },
      {
        type: 'code',
        code: `// Tiny RSA round-trip (toy values)\nconst { inverse_mod, power_mod } = Sage;\n\nconst p = 61n;\nconst q = 53n;\nconst n = p * q;\nconst phi = (p - 1n) * (q - 1n);\nconst e = 17n;\nconst d = inverse_mod(e, phi);\n\nconst message = 65n;\nconst cipher = power_mod(message, e, n);\nconst plain = power_mod(cipher, d, n);\n\nprint('n =', n);\nprint('cipher =', cipher);\nprint('plain =', plain);`,
      },
    ],
  },
  {
    id: 'ec',
    title: '7. Elliptic curves over finite fields',
    blocks: [
      {
        type: 'paragraph',
        text:
          'Elliptic curves power modern key exchange and signatures. The group law is defined geometrically but computed with finite field algebra.',
      },
      {
        type: 'paragraph',
        text:
          'We use the short Weierstrass form y^2 = x^3 + ax + b over a prime field GF(p). Points form a finite abelian group.',
      },
      {
        type: 'code',
        code: `// Elliptic curve taste\nconst { EllipticCurve, GF } = Sage;\n\nconst F = GF(97n);\nconst E = EllipticCurve(F, [2n, 3n]);\nconst P = E.point([3n, 6n]);\nconst Q = P.mul(20n);\n\nprint('P =', P.toString());\nprint('20P =', Q.toString());`,
      },
    ],
  },
];
