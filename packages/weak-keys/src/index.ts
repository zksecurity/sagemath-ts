/** Offline research screens, not a production cryptography or primality-proof API. */
import { Z_isanypower } from '../../parigp-ts/src/ifactor.js';
import { is_prime } from '../../sagemath-ts/src/arith/misc.js';
import {
  findSmallEcKeys,
  validateEcKeys,
} from '../../zksecurity-cheatsheets/src/attacks/ec-keys.js';
import { fermatFactor } from '../../zksecurity-cheatsheets/src/attacks/rsa-keys.js';
import { p256, p384, secp256k1 } from '../../zksecurity-cheatsheets/src/curves/params.js';

export const CHECKS = [
  'rsa-prime-modulus',
  'rsa-perfect-power',
  'rsa-fermat',
  'ec-small-scalar',
] as const;
export type Check = (typeof CHECKS)[number];
export type Status =
  | 'matched'
  | 'no-match'
  | 'not-run'
  | 'budget-exhausted'
  | 'malformed'
  | 'unsupported';
export type Confidence = 'exact' | 'probable' | 'bounded' | 'none';
export interface ScreenOptions {
  /** Default: prime-modulus and perfect-power only. Expensive screens are opt-in. */
  checks?: readonly Check[];
  /** Number of records eligible for work, including invalid records; default 32, hard maximum 256. */
  maxKeys?: number;
  /** Default 4096, hard maximum 8192. Larger inputs are unsupported, not silently truncated. */
  maxModulusBits?: number;
  /** Steps per RSA key; default 256, hard maximum 4096. Zero means no search. */
  fermatIterations?: number;
  /** Total Fermat steps allowed across the batch; default 8192, hard maximum 65536. */
  maxFermatSteps?: number;
  /** Look for +/-scalars through 2^bits inclusive; default 8, hard maximum 16. */
  ecTableBits?: number;
  /** Number of EC keys eligible for work; default 8, hard maximum 32. */
  maxEcKeys?: number;
}
export interface CheckResult {
  check: Check;
  status: Status;
  confidence: Confidence;
  /** Present only after work: declared search range, not recovered secret metadata. */
  scope?: { iterations: number } | { maxScalar: number };
}
export interface ScreenReport {
  schemaVersion: 1;
  records: Array<{ index: number; checks: CheckResult[] }>;
}

const CURVES = { p256, p384, secp256k1 };
function bounded(value: number | undefined, fallback: number, maximum: number): number {
  const n = value ?? fallback;
  if (!Number.isSafeInteger(n) || n < 0 || n > maximum)
    throw new RangeError('Invalid screening budget');
  return n;
}
function object(value: unknown): value is Record<string, unknown> {
  return typeof value === 'object' && value !== null && !Array.isArray(value);
}
function integer(value: unknown): bigint | null {
  if (typeof value === 'bigint') return value;
  // At most 8192 bits in hex; decimal strings are bounded too, before constructing a BigInt.
  if (
    typeof value !== 'string' ||
    value.length > 2470 ||
    !/^(?:[0-9]+|0x[0-9a-fA-F]+)$/.test(value)
  )
    return null;
  return BigInt(value);
}

/**
 * Screen explicit local public-key records. RSA: {type:'rsa',n,e} (bigint or integer strings).
 * EC: {type:'ec',curve:'p256'|'p384'|'secp256k1',point:SEC1Hex}.
 *
 * No input values, key identifiers, recovered factors/scalars, exception strings or upstream
 * recovery messages are forwarded. Results use input indices and fixed enums only. Callers
 * must still decide whether aggregate findings are appropriate to publish.
 *
 * `no-match` means the declared check scope completed, NOT universal key safety. A probable
 * prime is not certified. CPU bounds are counts/bit sizes, not interruptible wall deadlines.
 * This is an application adapter, with no SageMath counterpart; underlying math is reused.
 */
export function screenPublicKeys(
  inputs: readonly unknown[],
  options: ScreenOptions = {}
): ScreenReport {
  if (!Array.isArray(inputs) || inputs.length > 10000)
    throw new RangeError('At most 10000 records per report');
  const selected = options.checks ?? ['rsa-prime-modulus', 'rsa-perfect-power'];
  if (!Array.isArray(selected) || selected.some((c) => !CHECKS.includes(c)))
    throw new RangeError('Unknown screen');
  const maxKeys = bounded(options.maxKeys, 32, 256);
  const maxBits = bounded(options.maxModulusBits, 4096, 8192);
  const stepsPerKey = bounded(options.fermatIterations, 256, 4096);
  let stepsLeft = bounded(options.maxFermatSteps, 8192, 65536);
  const tableBits = bounded(options.ecTableBits, 8, 16);
  if (tableBits < 1) throw new RangeError('EC table bits must be positive');
  const maxEcKeys = bounded(options.maxEcKeys, 8, 32);
  let ecUsed = 0;
  return {
    schemaVersion: 1,
    records: inputs.map((input, index) => {
      // Validate each in-budget record once; do not echo arbitrary source fields.
      const record = object(input) ? input : null;
      let n: bigint | null = null;
      let e: bigint | null = null;
      let rsaState: 'malformed' | 'unsupported' | null = null;
      if (index < maxKeys && record?.type === 'rsa') {
        n = integer(record.n);
        e = integer(record.e);
        if (n === null || e === null || n <= 1n || e < 3n || !(e & 1n) || e >= n)
          rsaState = 'malformed';
        else if (n.toString(2).length > maxBits || e.toString(2).length > maxBits)
          rsaState = 'unsupported';
      }
      let ecState: 'malformed' | 'unsupported' | 'budget-exhausted' | null = null;
      let ecMatched = false;
      if (index < maxKeys && record?.type === 'ec' && selected.includes('ec-small-scalar')) {
        if (typeof record.curve !== 'string' || !Object.hasOwn(CURVES, record.curve))
          ecState = 'unsupported';
        else if (
          typeof record.point !== 'string' ||
          record.point.length > 200 ||
          !/^(?:0x)?[0-9a-fA-F]+$/.test(record.point)
        )
          ecState = 'malformed';
        else if (ecUsed >= maxEcKeys) ecState = 'budget-exhausted';
        else {
          ecUsed++;
          const curve = CURVES[record.curve as keyof typeof CURVES];
          if (!validateEcKeys(curve, [record.point])[0]!.valid) ecState = 'malformed';
          else {
            const bound = 1n << BigInt(tableBits);
            // A cached upstream table can be larger than requested. Filter recovered values
            // to our declared interval so process history cannot expand the advertised scope.
            ecMatched = findSmallEcKeys(curve, [record.point], { tableBits }).some(
              (result) => result.scalar <= bound || result.scalar >= curve.order - bound
            );
          }
        }
      }
      return {
        index,
        checks: CHECKS.map((check): CheckResult => {
          const result = (status: Status, confidence: Confidence = 'none'): CheckResult => ({
            check,
            status,
            confidence,
          });
          if (!selected.includes(check)) return result('not-run');
          if (index >= maxKeys) return result('budget-exhausted');
          if (!record || typeof record.type !== 'string') return result('malformed');
          if (record.type !== 'rsa' && record.type !== 'ec') return result('unsupported');
          if (check === 'ec-small-scalar') {
            if (record.type !== 'ec') return result('unsupported');
            if (ecState) return result(ecState);
            return {
              ...result(ecMatched ? 'matched' : 'no-match', ecMatched ? 'exact' : 'bounded'),
              scope: { maxScalar: 2 ** tableBits },
            };
          }
          if (record.type !== 'rsa') return result('unsupported');
          if (rsaState) return result(rsaState);
          // Non-null follows RSA input validation above; no partial keys reach math helpers.
          const modulus = n!;
          if (check === 'rsa-prime-modulus') {
            const passed = is_prime(modulus);
            return result(passed ? 'matched' : 'no-match', passed ? 'probable' : 'exact');
          }
          if (check === 'rsa-perfect-power')
            return result(Z_isanypower(modulus)[0] > 1 ? 'matched' : 'no-match', 'exact');
          const steps = Math.min(stepsPerKey, stepsLeft);
          if (!steps) return result('budget-exhausted');
          stepsLeft -= steps;
          const factors = fermatFactor(modulus, { maxIterations: steps });
          const matched =
            factors !== null &&
            factors[0] > 1n &&
            factors[1] > 1n &&
            factors[0] * factors[1] === modulus;
          return {
            ...result(matched ? 'matched' : 'no-match', matched ? 'exact' : 'bounded'),
            scope: { iterations: steps },
          };
        }),
      };
    }),
  };
}
