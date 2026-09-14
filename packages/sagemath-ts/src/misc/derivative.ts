import { OverflowError, ValueError } from '../errors.js';
import { Integer } from '../rings/integer_ring.js';

/**
 * Expand Sage's variable/count derivative argument sequence.
 * @see Deviation: Polynomial Derivative Protocol and Native Coefficients
 */
export function derivative_parse(args: readonly unknown[]): unknown[] {
  if (!args.length) return [null];
  if (args.length === 1 && Array.isArray(args[0])) return args[0];
  const output: unknown[] = [];
  let gotVariable = false,
    variable: unknown;
  for (const arg of args) {
    if (
      arg instanceof Integer ||
      typeof arg === 'bigint' ||
      typeof arg === 'boolean' ||
      (typeof arg === 'number' && Number.isInteger(arg))
    ) {
      const n = arg instanceof Integer ? arg.value : BigInt(arg);
      if (n < -(1n << 63n) || n >= 1n << 63n)
        throw new OverflowError('Python int too large to convert to C long');
      if (n < -(1n << 31n) || n >= 1n << 31n)
        throw new OverflowError('value too large to convert to int');
      if (n < 0n) throw new ValueError('derivative counts must be nonnegative');
      if (!gotVariable) variable = null;
      for (let i = 0; i < Number(n); i++) output.push(variable);
      gotVariable = false;
    } else {
      if (gotVariable) output.push(variable);
      gotVariable = true;
      variable = arg;
    }
  }
  if (gotVariable) output.push(variable);
  return output;
}
/**
 * Apply _derivative in the expanded variable order, retaining zero-count identity.
 * @see Deviation: Polynomial Derivative Protocol and Native Coefficients
 */
export function multi_derivative<F extends { _derivative(variable?: unknown): F }>(
  F: F,
  args: readonly unknown[]
): F {
  if (!args.length) return F._derivative();
  for (const arg of derivative_parse(args)) F = F._derivative(arg);
  return F;
}
