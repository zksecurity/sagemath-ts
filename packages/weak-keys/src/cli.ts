#!/usr/bin/env bun
/** Local JSON input only; no network, discovery, target connections or private-key output. */
import { closeSync, openSync, readSync } from 'node:fs';
import { CHECKS, type Check, type ScreenOptions, screenPublicKeys } from './index.js';

const HELP = `weak-keys --input FILE|- [--checks ${CHECKS.join(',')}]
  --max-keys N             eligible records (default 32; max 256)
  --max-modulus-bits N     RSA input bound (default 4096; max 8192)
  --fermat-iterations N    steps per RSA key (default 256; max 4096)
  --max-fermat-steps N     total steps (default 8192; max 65536)
  --ec-table-bits N        +/-small-scalar bound 2^N (default 8; max 16)
  --max-ec-keys N          EC record budget (default 8; max 32)
Input: JSON array of {type:'rsa',n:'decimal/0xhex',e:'decimal/0xhex'} or
{type:'ec',curve:'p256'|'p384'|'secp256k1',point:'SEC1hex'}.
Default checks: rsa-prime-modulus,rsa-perfect-power. Output: allowlisted JSON report.
Input limit 1 MiB. No network. No-match does not imply a safe key.
Timestamp profiles: see packages/weak-keys/timestampseeds (Python CLI).`;

function readInput(path: string): unknown {
  const fd = path === '-' ? 0 : openSync(path, 'r');
  try {
    const chunks: Buffer[] = [];
    let total = 0;
    for (;;) {
      const chunk = Buffer.alloc(65536);
      const size = readSync(fd, chunk, 0, chunk.length, null);
      if (!size) break;
      total += size;
      if (total > 1048576) throw new RangeError('Input limit');
      chunks.push(chunk.subarray(0, size));
    }
    return JSON.parse(Buffer.concat(chunks).toString('utf8'));
  } finally {
    if (path !== '-') closeSync(fd);
  }
}

/** Main returns an exit status and never echoes untrusted records or error messages. */
export function main(args: readonly string[]): number {
  if (!args.length || args.includes('--help')) {
    console.log(HELP);
    return 0;
  }
  try {
    let input: string | undefined;
    const options: ScreenOptions = {};
    const numeric = {
      '--max-keys': 'maxKeys',
      '--max-modulus-bits': 'maxModulusBits',
      '--fermat-iterations': 'fermatIterations',
      '--max-fermat-steps': 'maxFermatSteps',
      '--ec-table-bits': 'ecTableBits',
      '--max-ec-keys': 'maxEcKeys',
    } as const;
    const seen = new Set<string>();
    for (let i = 0; i < args.length; i += 2) {
      const flag = args[i]!;
      const value = args[i + 1];
      if (value === undefined || seen.has(flag)) throw new Error('Option');
      seen.add(flag);
      if (flag === '--input') input = value;
      else if (flag === '--checks') {
        const checks = value.split(',');
        if (checks.some((check) => !CHECKS.includes(check as Check))) throw new Error('Check');
        options.checks = checks as Check[];
      } else if (Object.hasOwn(numeric, flag) && /^\d+$/.test(value)) {
        options[numeric[flag as keyof typeof numeric]] = Number(value);
      } else throw new Error('Option');
    }
    if (!input) throw new Error('Input');
    const keys = readInput(input);
    if (!Array.isArray(keys)) throw new Error('Array required');
    console.log(JSON.stringify(screenPublicKeys(keys, options)));
    return 0;
  } catch {
    console.error('weak-keys: invalid input or options; use --help');
    return 2;
  }
}

if (import.meta.main) process.exitCode = main(process.argv.slice(2));
