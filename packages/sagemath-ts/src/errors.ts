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
