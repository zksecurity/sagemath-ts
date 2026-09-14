import {
  sha256,
  hmac_sha256,
  DeriveKey,
  RandomStream,
  salsa20_core,
  salsa20_apply,
  salsa20_init,
} from '../../../packages/ntl-ts/src/ZZ.js';
export function ntl_random_stream(
  op: bigint,
  n: bigint,
  av: bigint[],
  bv: bigint[],
  commands: bigint[]
): string {
  const a = Uint8Array.from(av, Number),
    b = Uint8Array.from(bv, Number);
  const hex = (data: Uint8Array) =>
    Array.from(data, (x) => x.toString(16).padStart(2, '0')).join('');
  const words = (state: Uint32Array) => Array.from(state, String);
  let out: unknown;
  if (op === 0n) out = hex(sha256(a, Number(n)));
  else if (op === 1n) out = hex(hmac_sha256(b, a, Number(n)));
  else if (op === 2n) out = hex(DeriveKey(a, Number(n)));
  else if (op === 7n) {
    // A view with a reported native pointer length; only its low 32 bits are read.
    Object.defineProperty(a, 'length', { value: Number(commands[0]) });
    out = hex(hmac_sha256(b, a, Number(n)));
  } else if (op === 3n) {
    const streams = [new RandomStream(a), new RandomStream(a), new RandomStream(a)],
      trace: unknown[] = [];
    for (let i = 0; i < commands.length; i += 3) {
      const action = commands[i],
        slot = Number(commands[i + 1]),
        arg = commands[i + 2]!;
      try {
        if (action === 0n) trace.push([null, hex(streams[slot]!.get(Number(arg)))]);
        else {
          if (action === 1n) streams[slot]!.set_nonce(arg);
          else if (action === 2n) streams[slot] = new RandomStream(streams[Number(arg)]!);
          else if (action === 3n) streams[slot] = streams[Number(arg)]!;
          else if (action === 4n) streams[slot]!.assign(streams[Number(arg)]!);
          else throw new Error('unknown stream command');
          trace.push([null, null]);
        }
      } catch (e) {
        trace.push([(e as Error).message, null]);
      }
    }
    out = trace;
  } else if (op === 4n || op === 5n) {
    const state = Uint32Array.from(av, Number);
    if (op === 4n) {
      salsa20_core(state);
      out = words(state);
    } else {
      const block = salsa20_apply(state);
      out = [words(state), words(block)];
    }
  } else if (op === 6n) out = words(salsa20_init(a));
  else throw new Error('unknown NTL random stream operation');
  return JSON.stringify(out);
}
