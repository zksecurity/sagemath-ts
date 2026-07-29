/**
 * @module sage/errors
 * @description Custom error classes matching SageMath's error types
 */

export class ValueError extends Error {
  override name = 'ValueError';

  constructor(message: string) {
    super(message);
  }
}

export const TypeError = globalThis.TypeError;

export class ZeroDivisionError extends Error {
  override name = 'ZeroDivisionError';

  constructor(message: string = 'division by zero') {
    super(message);
  }
}

export class NotImplementedError extends Error {
  override name = 'NotImplementedError';

  constructor(message: string = 'not implemented') {
    super(message);
  }
}

export class ArithmeticError extends Error {
  override name = 'ArithmeticError';

  constructor(message: string) {
    super(message);
  }
}

export class PrecisionError extends Error {
  override name = 'PrecisionError';

  constructor(message: string = 'precision error') {
    super(message);
  }
}

/**
 * Python's built-in `RuntimeError`, raised by SageMath where no more specific
 * exception applies (e.g. `discrete_gaussian_lattice.py:570`).
 */
export class RuntimeError extends Error {
  override name = 'RuntimeError';

  constructor(message: string) {
    super(message);
  }
}

/**
 * Python's built-in `IndexError`.
 *
 * Raised by SageMath when a coefficient is asked for beyond what is known,
 * e.g. `power_series_ring_element.pyx` / `laurent_series_ring_element.pyx`
 * `__getitem__` past the series precision.
 */
export class IndexError extends Error {
  override name = 'IndexError';

  constructor(message: string) {
    super(message);
  }
}

/**
 * Python's built-in `AssertionError`, raised by a failing `assert` statement.
 */
export class AssertionError extends Error {
  override name = 'AssertionError';

  constructor(message: string = '') {
    super(message);
  }
}

/**
 * cypari2's `PariError`, i.e. an error raised inside PARI itself.
 *
 * SageMath lets these propagate unchanged, so the message must be PARI's own
 * (e.g. `Qfb: sorry, negative definite t_QFB is not yet implemented`).
 */
export class PariError extends Error {
  override name = 'PariError';

  constructor(message: string) {
    super(message);
  }
}
