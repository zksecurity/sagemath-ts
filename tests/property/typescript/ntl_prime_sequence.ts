import { PrimeSeq } from '../../../packages/ntl-ts/src/ZZ.js';
export function ntl_prime_sequence(commands: bigint[]): string {
  const seq = [new PrimeSeq(), new PrimeSeq()],
    out: unknown[] = [];
  for (let i = 0; i < commands.length; i += 3) {
    const op = commands[i],
      slot = Number(commands[i + 1]),
      arg = commands[i + 2]!;
    if (op === 0n) out.push(Array.from({ length: Number(arg) }, () => seq[slot]!.next()));
    else if (op === 1n) {
      seq[slot]!.reset(arg);
      out.push(null);
    } else if (op === 2n) {
      seq[slot] = new PrimeSeq();
      out.push(null);
    } else throw new Error('unknown NTL prime sequence operation');
  }
  return JSON.stringify(out, (_, v) => (typeof v === 'bigint' ? String(v) : v));
}
