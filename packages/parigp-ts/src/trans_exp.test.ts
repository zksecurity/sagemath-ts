import { expect, test } from 'bun:test';
import { exp1r_abs } from './trans1.js';
import { PariError } from './errors.js';

test('exp1r_abs preserves every native exponent bit during repeated doubling', () => {
  // Exact frames from the bundled PARI 2.18.1 exp1r_abs, also exercised by the
  // shared Python/TypeScript pari_real_exp1 cases at multiple precisions.
  for (const s of [-1, 1] as const) {
    expect(exp1r_abs({s,e:54,p:384,m:19701003098197239606139520050071806902539869635232723333974146702122860885748605305707133127442457820403313995153408n})).toEqual({
      s:1,e:25989283394227192n,p:384,m:25402264293164475225123981458011408035153460912115219070077427635023958917720485126445057062292716840191275494591202n,
    });
    expect(exp1r_abs({s,e:59,p:512,m:1n << 511n})).toEqual({
      s:1,e:831657068615270155n,p:512,m:11150356226034666957404441079638023940694805724257669351967712767525258533424313847170998844709680992299767476244436355943581751173069022031140943750582060n,
    });
  }
});

test('exp1r_abs reports native exponent overflow and preserves tiny inputs', () => {
  expect(() => exp1r_abs({s:1,e:61,p:512,m:1n << 511n}))
    .toThrow(new PariError('overflow in expo()'));
  expect(exp1r_abs({s:-1,e:-64,p:64,m:1n << 63n}))
    .toEqual({s:1,e:-64n,p:64,m:1n << 63n});
});
