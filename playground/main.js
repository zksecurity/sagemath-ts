import { createEditor } from './editor.js';
import { createRunner, resolveSage } from './runner.js';

const textarea = document.getElementById('code');
const output = document.getElementById('output');
const status = document.getElementById('status');
const runtime = document.getElementById('runtime');
const runBtn = document.getElementById('run');
const resetBtn = document.getElementById('reset');
const copyBtn = document.getElementById('copy');

const examples = {
  gcd: `// GCD + Euclid\nconst { gcd } = Sage;\n\nconst pairs = [\n  [48n, 18n],\n  [99n, 78n],\n  [252n, 198n],\n];\n\nfor (const [a, b] of pairs) {\n  print('gcd(' + a + ', ' + b + ') =', gcd(a, b));\n}`,
  modinv: `// Modular inverse demo\nconst { inverse_mod, power_mod } = Sage;\n\nconst a = 37n;\nconst n = 101n;\nconst inv = inverse_mod(a, n);\n\nprint('a =', a, 'mod', n);\nprint('inverse =', inv);\nprint('check:', power_mod(a, 1n, n), '*', inv, '% n =', (a * inv) % n);`,
  rsa: `// Tiny RSA walkthrough (toy values)\nconst { inverse_mod, power_mod } = Sage;\n\nconst p = 61n;\nconst q = 53n;\nconst n = p * q;\nconst phi = (p - 1n) * (q - 1n);\nconst e = 17n;\nconst d = inverse_mod(e, phi);\n\nconst message = 65n;\nconst cipher = power_mod(message, e, n);\nconst plain = power_mod(cipher, d, n);\n\nprint('n =', n);\nprint('e =', e);\nprint('d =', d);\nprint('cipher =', cipher);\nprint('plain =', plain);`,
  ec: `// Elliptic curve taste\nconst { EllipticCurve, GF } = Sage;\n\nconst F = GF(97n);\nconst E = EllipticCurve(F, [2n, 3n]);\nconst P = E.point([3n, 6n]);\nconst Q = P.mul(20n);\n\nprint('P =', P.toString());\nprint('20P =', Q.toString());`,
};

const editor = createEditor(textarea);
const runner = createRunner({ outputEl: output, statusEl: status, runtimeEl: runtime });

function loadExample(name) {
  editor.setValue(examples[name] || examples.gcd);
  status.textContent = 'Idle';
}

runBtn.addEventListener('click', () => runner.run(editor.getValue()));
resetBtn.addEventListener('click', () => loadExample('gcd'));
copyBtn.addEventListener('click', async () => {
  await navigator.clipboard.writeText(editor.getValue());
  copyBtn.textContent = 'Copied';
  setTimeout(() => (copyBtn.textContent = 'Copy'), 1200);
});

document.querySelectorAll('[data-example]').forEach((button) => {
  button.addEventListener('click', () => loadExample(button.dataset.example));
});

document.addEventListener('keydown', (event) => {
  if ((event.metaKey || event.ctrlKey) && event.key === 'Enter') {
    runner.run(editor.getValue());
  }
});

loadExample('gcd');
resolveSage(runtime);
