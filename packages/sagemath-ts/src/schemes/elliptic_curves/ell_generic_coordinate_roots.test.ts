import { nativeFixtures as loadLiveNative } from "../../../../../tests/property/native-live.mjs";


import { expect, test } from 'bun:test';
import { execFileSync } from 'node:child_process';
import { GF } from '../../rings/finite_rings/finite_field_constructor.js';
import { GFpn } from '../../rings/finite_rings/finite_field_extension.js';
import { QQ } from '../../rings/rational_field.js';
import { Polynomial as IsomorphismTracePolynomial } from '../../rings/polynomial/polynomial_element.js';
import { EllipticCurve } from './constructor.js';
import { getrand as coordinateGetrand, setrand as coordinateSetrand } from '@sagemath-ts/parigp-ts';
const native = (await loadLiveNative(import.meta.url, "./ell_generic_coordinate_roots.native.json.gz")) as { args: (string)[]; error: null; errorType: null; function: string; result: string; seed: number }[];
type Any = any;
const functions: Record<string, (...args: Any[]) => string> = {};
const field = (p: bigint): Any => (p === 0n ? QQ : GF(p));
functions.ec_coordinate_roots = (
  p: bigint,
  degree: bigint,
  modulus: bigint[],
  coefficients: Any[],
  coordinate: Any,
  operation: bigint,
  seed: bigint
) => {
  const trace: Any[] = [];
  let result: Any;
  coordinateSetrand(seed);
  try {
    const K: Any = degree > 1n ? GFpn(p, Number(degree), modulus.slice(0, -1) as Any, 'a') : field(p);
    const decode = (v: Any) =>
      degree > 1n
        ? K.fromInteger(BigInt(v))
        : Array.isArray(v)
          ? K.__call__(v[0]).div(K.__call__(v[1]))
          : K.__call__(v);
    const encode = (v: Any) => (degree > 1n ? String(v.integer_representation()) : String(v));
    const E = EllipticCurve(K, coefficients.map(decode) as Any);
    const proto: Any = Object.getPrototypeOf(K.zero()),
      poly: Any = IsomorphismTracePolynomial.prototype;
    const square = proto.is_square,
      sqrt = proto.sqrt,
      roots = poly.roots;
    let depth = 0;
    try {
      proto.is_square = function (this: Any) {
        if (!depth) trace.push(['is_square', encode(this)]);
        depth++;
        try {
          return square.call(this);
        } finally {
          depth--;
        }
      };
      proto.sqrt = function (this: Any, options: Any) {
        if (!depth) trace.push(['sqrt', encode(this), options?.all ?? false]);
        depth++;
        try {
          return sqrt.call(this, options);
        } finally {
          depth--;
        }
      };
      poly.roots = function (this: Any, options: Any) {
        if (!depth) trace.push(['roots', this.coeffs.map(encode), options?.multiplicities ?? true]);
        depth++;
        try {
          return roots.call(this, options);
        } finally {
          depth--;
        }
      };
      let value: Any;
      if (operation === 0n) value = E.is_x_coord(decode(coordinate));
      else if (operation === 4n) value = String(E);
      else if (operation === 3n) value = (E.montgomery_model() as Any).a_invariants().map(encode);
      else {
        let pts: Any =
          operation === 1n ? E.lift_x(decode(coordinate), true) : [E.lift_x(decode(coordinate))];
        value = pts.map((P: Any) => [
          encode(P.x()),
          encode(P.y()),
          P.curve === E,
          P.xyz().every((v: Any) => p === 0n || v.parent === K),
        ]);
      }
      result = { value };
    } finally {
      proto.is_square = square;
      proto.sqrt = sqrt;
      poly.roots = roots;
    }
  } catch (e) {
    result = { error: (e as Error).name, message: (e as Error).message };
  }
  result.calls = trace;
  if (degree > 1n && p !== 2n && operation < 3n) result.state = String(coordinateGetrand());
  return JSON.stringify(result);
};

const watched = new Set<number>([29800000]);
for (const p of [
  '7',
  '3',
  '2',
  '2305843009213693951',
  '18446744073709551557',
  '618970019642690137449562111',
]) {
  const row = native.find(
    (r) =>
      r.seed >= 29600000 &&
      r.seed < 29700000 &&
      r.args[0] === p &&
      r.args[5] === '1' &&
      (p !== '2' || r.args[1] === '130') &&
      JSON.parse(r.result).value?.length === 2
  );
  if (row) watched.add(row.seed);
}
// Binary fields can have a single root; retain a degree-130 watch explicitly.
watched.add(
  native.find(
    (r) =>
      r.seed >= 29600000 &&
      r.seed < 29700000 &&
      r.args[0] === '2' &&
      r.args[1] === '130' &&
      r.args[5] === '1'
  )!.seed
);
for (const row of native) {
  test(`native curve coordinate roots ${row.seed}`, () => {
    const raw = row.args.map((arg) =>
      arg.startsWith('[') ? JSON.parse(arg.replace(/-?\d+/g, '"$&"')) : arg
    );
    const decode = (v: Any): Any => (Array.isArray(v) ? v.map(decode) : BigInt(v));
    if (watched.has(row.seed)) {
      const script = `
    import {GF} from '${import.meta.dir}/../../rings/finite_rings/finite_field_constructor.ts';
    import {GFpn} from '${import.meta.dir}/../../rings/finite_rings/finite_field_extension.ts';
    import {QQ} from '${import.meta.dir}/../../rings/rational_field.ts';
    import {Polynomial as IsomorphismTracePolynomial} from '${import.meta.dir}/../../rings/polynomial/polynomial_element.ts';
    import {EllipticCurve} from '${import.meta.dir}/constructor.ts';
    import {getrand as coordinateGetrand,setrand as coordinateSetrand} from '@sagemath-ts/parigp-ts';
    const field=p=>p===0n?QQ:GF(p);
    const decode=v=>Array.isArray(v)?v.map(decode):BigInt(v);
    console.log((${functions.ec_coordinate_roots!.toString()})(...${JSON.stringify(raw)}.map(decode)));`;
      expect(
        execFileSync(process.execPath, ['--eval', script], {
          encoding: 'utf8',
          timeout: 10000,
        }).trim()
      ).toBe(row.result);
    } else expect(functions.ec_coordinate_roots!(...raw.map(decode))).toBe(row.result);
  }, 15000);
}
