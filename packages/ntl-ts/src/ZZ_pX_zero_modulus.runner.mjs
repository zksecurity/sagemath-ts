// Run the potentially nonterminating regression in a process the unit test can stop.
import { ZZ_pXModulus, build, rem } from './ZZ_pX.ts';
import { FFTPrimeContext } from './FFT.ts';
import { RandomStream } from './ZZ.ts';
const parse = (s) => s.slice(1, -1).trim() ? s.slice(1, -1).split(',').map((x) => BigInt(x.trim())) : [];
const hex = (a) => Array.from(a, (x) => x.toString(16).padStart(2, '0')).join('');
const bits = (v) => {
  const view = new DataView(new ArrayBuffer(8));
  view.setFloat64(0, v, false);
  return view.getBigUint64(0, false).toString(16).padStart(16, '0');
};
const info = (v) => [v.q, bits(v.qrecip), v.RootTable[0], v.RootTable[1], v.TwoInvTable];
const stringify = (v) => JSON.stringify(v, (_, x) => typeof x === 'bigint' ? String(x) : x);
const records = JSON.parse(await new Response(Bun.stdin.stream()).text());
const results = records.map((args) => {
  let result = null, error = null, errorType = null;
  try {
    const [key, p, lengths, coefficients, commands] = [parse(args[0]), BigInt(args[1]), parse(args[2]), parse(args[3]), parse(args[4])];
    const context = new FFTPrimeContext(), stream = new RandomStream(Uint8Array.from(key, Number));
    const F = new ZZ_pXModulus(null, p, {context,stream});
    let offset = 0;
    const polys = lengths.map((n) => {const a = coefficients.slice(offset, offset + Number(n)); offset += Number(n); return a;});
    const trace = [];
    for (let i=0; i<commands.length; i+=2) {
      let value = null, kind = null, message = null;
      try {
        const op = Number(commands[i]), a = polys[Number(commands[i+1])];
        if (op === 0) build(F, a);
        else if (op === 1) value = rem(a, F);
        else throw new Error('unexpected constant-modulus regression operation');
      } catch(e) {kind=e.name;message=e.message;}
      trace.push([kind,message,value,BigInt(F.n),F.f.slice(),BigInt(F.UseFFT),Array.from({length:context.length()},(_,j)=>info(context.get(j))),hex(new RandomStream(stream).get(64))]);
    }
    result = stringify(trace);
  } catch(e) {error=e.message;errorType=e.name;}
  return {result,error,errorType};
});
process.stdout.write(JSON.stringify(results));
