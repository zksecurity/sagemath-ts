import { describe, expect, spyOn, test } from 'bun:test';
import * as kernels from './_modular_sqrt.js';
import * as bb from './bb_group.js';
import * as ff from './ff.js';

describe('PARI modular square-root fidelity', () => {
  test('native zero/composite boundaries terminate with the original result or error', () => {
    expect(ff.Fp_sqrt(0n, 0n)).toBe(0n);
    expect(() => ff.Fp_sqrt(1n, 0n)).toThrow('impossible inverse in dvmdii: 0.');
    for (const p of [9n, 25n, 49n, 81n])
      expect(() => ff.Fp_sqrt(1n, p)).toThrow(`not a prime number in Fl_nonsquare: ${p}.`);
    expect(ff.Fp_sqrt(2n, -17n)).toBe(6n);
  });

  test('word dispatch and signed-small shortcuts retain their native routes', () => {
    const word = spyOn(kernels, 'wordSquareRoot');
    const power = spyOn(ff, 'Fp_pow');
    try {
      expect(ff.Fp_sqrt(1n, 17n)).toBe(1n);
      expect(word).toHaveBeenCalledTimes(1);
      expect(power).toHaveBeenCalledTimes(2);
      word.mockClear();
      power.mockClear();
      expect(ff.Fp_sqrt(1n, 36893488147419103397n)).toBe(1n);
      expect(word).toHaveBeenCalledTimes(0);
      expect(power).toHaveBeenCalledTimes(0);
    } finally {
      word.mockRestore();
      power.mockRestore();
    }
  });

  test('Cipolla uses native fused powering above its cost threshold', () => {
    const fold = spyOn(bb, 'gen_pow_fold');
    try {
      expect(ff.Fp_sqrt(4099276460830869277n, 36893488147477823489n)).toBe(12297829382492607829n);
      expect(fold).toHaveBeenCalledTimes(0);
      expect(ff.Fp_sqrt(8198552921741896363n, 36893488147838533633n)).toBe(12297829382612844544n);
      expect(fold).toHaveBeenCalledTimes(1);
      expect(fold.mock.calls[0]![1]).toBe(18446744073919266816n);
    } finally {
      fold.mockRestore();
    }
  });
});
