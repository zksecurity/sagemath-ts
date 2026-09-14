import { nativeFixtures as loadLiveNative } from "../../../../../tests/property/native-live.mjs";
import { execFileSync } from 'node:child_process';
import { expect, test } from 'bun:test';
const fixtures = await loadLiveNative(import.meta.url, "./division_points_order.native.json");
import { GF } from '../../rings/finite_rings/finite_field_constructor.js';
import { EllipticCurve } from './constructor.js';
import { division_points } from './ell_point.js';

for (const row of fixtures) {
  test(`bundled Sage division points ${row.seed}`, () => {
    const parse = (s: string) =>
      s.slice(1, -1).trim()
        ? s
            .slice(1, -1)
            .split(',')
            .map((c) => BigInt(c.trim()))
        : [];
    const K = GF(BigInt(row.args[0]!)),
      E = EllipticCurve(K, parse(row.args[1]!) as [bigint, bigint, bigint, bigint, bigint]);
    const point = parse(row.args[2]!),
      m = BigInt(row.args[3]!);
    const P = point.length ? E.point([K.__call__(point[0]!), K.__call__(point[1]!)]) : E.zero();
    const before = P.is_zero() ? [] : [String(P.x()), String(P.y())];
    const evaluate = () => {
      const value = division_points(P, m).map((Q) =>
        Q.is_zero()
          ? [['0', '1', '0'], true]
          : [[String(Q.x()), String(Q.y()), '1'], Q.x().parent === K && Q.y().parent === K]
      );
      return JSON.stringify({ value });
    };
    let result: string;
    if (m >= 100n) {
      // This degree-5100 regression must reduce to gcd(f,x^q-x), as Sage does.
      // A timeout is failure, never an expected comparison result.
      const program = `import {GF} from ${JSON.stringify(new URL('../../rings/finite_rings/finite_field_constructor.ts', import.meta.url).pathname)};
        import {EllipticCurve} from ${JSON.stringify(new URL('./constructor.ts', import.meta.url).pathname)};
        import {division_points} from ${JSON.stringify(new URL('./ell_point.ts', import.meta.url).pathname)};
        const K=GF(${row.args[0]}n),E=EllipticCurve(K,${JSON.stringify(parse(row.args[1]!).map(String))}.map(BigInt));
        const m=${m}n,point=${JSON.stringify(point.map(String))}.map(BigInt);
        const P=point.length?E.point([K.__call__(point[0]),K.__call__(point[1])]):E.zero();
        console.log((${evaluate.toString()})());`;
      result = execFileSync(process.execPath, ['-e', program], {
        encoding: 'utf8',
        timeout: 10000,
      }).trim();
    } else result = evaluate();
    expect({ result, error: null, errorType: null }).toEqual({
      result: row.result,
      error: row.error,
      errorType: row.errorType,
    });
    expect(P.is_zero() ? [] : [String(P.x()), String(P.y())]).toEqual(before);
  }, 30000);
}
