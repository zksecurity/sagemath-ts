import * as Sage from 'sagemath-ts';

// Expose a global for the playground runtime.
if (typeof globalThis !== 'undefined') {
  (globalThis as typeof globalThis & { Sage?: typeof Sage }).Sage = Sage;
}

export {};
