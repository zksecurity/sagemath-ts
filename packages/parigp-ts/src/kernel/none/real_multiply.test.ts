import { expect, test } from 'bun:test';
import { mulir, mulri, mulrr, sqrr, type MpReal } from '../../qfb.js';

const boundary = (p: number): MpReal => ({ s: 1, e: 0, p, m: (5n << BigInt(p - 3)) - 3n });

test('native small square keeps the truncated guard sum', () => {
  expect(sqrr(boundary(192))).toEqual({ s: 1, e: 0, p: 192, m: (25n << 187n) - 8n });
});

test('native pointer identity chooses the square threshold', () => {
  const x = boundary(832);
  expect(mulrr(x, x)).toEqual({ s: 1, e: 0, p: 832, m: (25n << 827n) - 7n });
  expect(mulrr(x, { ...x })).toEqual({ s: 1, e: 0, p: 832, m: (25n << 827n) - 8n });
});

test('native mixed integer products use the same truncated kernel', () => {
  const n = 4638679697461769301717983489047437979229256876716044252275n;
  const y: MpReal = { s: 1, e: 0, p: 192, m: 3297503506553735452683461125856361419140298496881443840374n };
  const expected: MpReal = { s: 1, e: 191, p: 192, m: 4873606709902253955708965301927856626279744602156248068738n };
  expect(mulir(n, y)).toEqual(expected);
  expect(mulri(y, n)).toEqual(expected);
});

test('word product rounding carries into the exponent', () => {
  const y: MpReal = { s: 1, e: 0, p: 64, m: ((1n << 66n) - 1n) / 7n };
  expect(mulir(7n, y)).toEqual({ s: 1, e: 3, p: 64, m: 1n << 63n });
});
