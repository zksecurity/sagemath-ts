import { nativeFixtures as loadLiveNative } from "../../../../../tests/property/native-live.mjs";
import { expect, test } from 'bun:test';
import { QuaternionAlgebra } from './quaternion_algebra.js';
import { Rational } from '../../rings/rational.js';
const fixtures = await loadLiveNative(import.meta.url, "./lll_reduction.native.json");
import resources from './lll_resources.native.json' with { type: 'json' };

function operation(name: string, args: string[]): string {
  if (name === 'quat_order_isomorphism') {
    const scale = BigInt(args[0]!),
      A = QuaternionAlgebra(-1n, -19n),
      [i, j, k] = A.gens();
    const O0 = A.quaternion_order([
      A.one(),
      i,
      i.add(j).scalar_mul('1/2'),
      A.one().add(k).scalar_mul('1/2'),
    ]);
    let O1 = A.quaternion_order([
      A.one(),
      i.scalar_mul(667n),
      A.__call__('1/2').add(j.scalar_mul('1/2')).add(i.scalar_mul(9n)),
      i.scalar_mul('222075/2').add(j.scalar_mul(333n)).add(k.scalar_mul('1/2')).scalar_mul('1/667'),
    ]);
    if (scale !== 0n) {
      const a = A.__call__([1n, scale, scale + 1n, 2n * scale - 1n]);
      O1 = A.quaternion_order(O1.basis().map((x) => a.inverse().mul(x).mul(a)));
    }
    const gamma = O0.isomorphism_to(O1, { conjugator: true }) as ReturnType<typeof A.__call__>;
    return [gamma, ...A.gens().map((x) => gamma.inverse().mul(x).mul(gamma))]
      .map((x) => x.list().map(String).join(','))
      .join(';');
  }
  const a = BigInt(args[0]!),
    b = BigInt(args[1]!),
    den = BigInt(args[2]!),
    op = BigInt(args[4]!);
  const flat = args[3]!
    .slice(1, -1)
    .split(',')
    .map((x) => BigInt(x.trim()));
  const A = QuaternionAlgebra(a, b),
    I = A.ideal(
      Array.from({ length: 4 }, (_, i) =>
        A.__call__(flat.slice(i * 4, i * 4 + 4)).scalar_mul(new Rational(1n, den))
      )
    );
  if (op === 0n)
    return I.reduced_basis()
      .map((x) => x.list().map(String).join(','))
      .join(';');
  if (op === 1n) return I.minimal_element().list().map(String).join(',');
  if (op === 2n)
    return Array.from({ length: 4 }, (_, i) =>
      Array.from({ length: 4 }, (_, j) => String(I.quadratic_form().get(i, j)))
    )
      .flat()
      .join(',');
  return I.theta_series_vector(12).join(',');
}
for (const [i, f] of fixtures.entries())
  test(`original quaternion basis, form and conjugator ${i}`, () => {
    let result: string | null = null,
      error: string | null = null,
      errorType: string | null = null;
    try {
      result = operation(f.function, f.args);
    } catch (e) {
      error = (e as Error).message;
      errorType = (e as Error).name;
    }
    expect({ result, error, errorType }).toEqual({
      result: f.result,
      error: f.error,
      errorType: f.errorType,
    });
  });
// Native PARI also fails on these forms (stack exhaustion); these are explicitly
// resource-boundary tests, not equal-error or successful comparative assertions.
for (const [i, f] of resources.entries())
  test(`native and host resource boundary for indefinite Gram ${i}`, () => {
    expect(f.native.error).toContain('PARI stack overflows');
    let error: Error | undefined;
    try {
      operation(f.port.function, f.port.args);
    } catch (e) {
      error = e as Error;
    }
    expect(error?.name).toBe(f.port.errorType);
    expect(error?.message).toBe(f.port.error);
  });
