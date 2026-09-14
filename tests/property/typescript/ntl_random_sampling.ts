import {
  RandomStream,
  RandomWord,
  VectorRandomWord,
  RandomBits_long,
  RandomBits_ulong,
  RandomLen_long,
  RandomBits,
  RandomBits_ZZ,
  RandomLen,
  RandomLen_ZZ,
  RandomBnd,
} from '../../../packages/ntl-ts/src/ZZ.js';
export function ntl_random_sampling(key: bigint[], cmd: bigint[]): string {
  const a = Uint8Array.from(key, Number),
    s = [new RandomStream(a), new RandomStream(a), new RandomStream(a)],
    out: unknown[] = [];
  const hex = (b: Uint8Array) => Array.from(b, (x) => x.toString(16).padStart(2, '0')).join('');
  for (let i = 0; i < cmd.length; i += 3) {
    const op = cmd[i],
      slot = Number(cmd[i + 1]),
      arg = cmd[i + 2]!,
      stream = s[slot]!;
    let value: unknown = null;
    try {
      if (op === 0n) value = RandomWord(stream);
      else if (op === 1n) value = VectorRandomWord(Number(arg), stream);
      else if (op === 2n) value = RandomBits_long(Number(arg), stream);
      else if (op === 3n) value = RandomBits_ulong(Number(arg), stream);
      else if (op === 4n) value = RandomLen_long(Number(arg), stream);
      else if (op === 5n) value = RandomBits(Number(arg), stream);
      else if (op === 6n) value = RandomLen(Number(arg), stream);
      else if (op === 7n) value = RandomBnd(arg, stream);
      else if (op === 8n) value = RandomBits_ZZ(Number(arg), stream);
      else if (op === 9n) value = RandomLen_ZZ(Number(arg), stream);
      else if (op === 10n) value = RandomBnd(arg, stream, { word: true });
      else if (op === 11n) value = hex(stream.get(Number(arg)));
      else if (op === 12n) stream.set_nonce(arg);
      else if (op === 13n) s[slot] = new RandomStream(s[Number(arg)]!);
      else if (op === 14n) stream.assign(s[Number(arg)]!);
      else if (op === 15n) s[slot] = s[Number(arg)]!;
      else throw new Error('unknown integer sampling operation');
      out.push([null, null, value]);
    } catch (e) {
      out.push([(e as Error).name, (e as Error).message, null]);
    }
  }
  for (const stream of s) out.push([null, null, hex(stream.get(64))]);
  return JSON.stringify(out, (_, v) => (typeof v === 'bigint' ? String(v) : v));
}
