import { expect, test } from 'bun:test';
import { logr_abs, mplog2, abpq_init, abpq_sum } from './trans1.js';
import { atanhuu } from './trans2.js';
import { itor, real_0_bit, sqrtr as qfbSqrt } from './qfb.js';
import { sqrtr, expr, rtodbl } from './buch.js';
import { PariError } from './errors.js';

test('native logarithm increases precision during its series recurrence', () => {
  expect(logr_abs({s:1,e:-1,p:384,m:39402006196394479208099832594110562112495958950223682931689384367727332642107873841363075178211817256653219393729825n})).toEqual({s:-1,e:-64,m:38546385647358835747933002318264302861693239821253018900752036738380079559949051445849012752476470991646796827996179n,p:384});
  expect(logr_abs({s:1,e:-1,p:512,m:13407807929942597098315134169554226086712180169958699491849565324361474642630198936489781267914149421850622757523434859188431596408010985175227994749046169n})).toEqual({s:-1,e:-64,m:11611218466438289295898820582638552556854148904254881001502724220891367383338435628051014700696110719377759153233144914325614438583206121995213717364481417n,p:512});
});

test('Buchmann square roots retain native complex results and zero accuracy', () => {
  expect(sqrtr(itor(-4n, 64))).toEqual({ re: 0n, im: itor(2n, 64) });
  expect(sqrtr(real_0_bit(2 ** 40 + 1))).toEqual(real_0_bit(2 ** 39));
  const x = {s:1 as const,e:0,p:64,m:(1n << 63n) + 1n};
  expect(sqrtr(x)).toEqual({s:1,e:0,p:64,m:(1n << 63n) + 1n});
  expect(sqrtr(x)).toEqual(qfbSqrt(x));
});

test('Buchmann exponential retains native zero accuracy and BigInt exponents', () => {
  expect(expr(real_0_bit(1))).toEqual({s:0,e:1n,p:0,m:0n});
  expect(expr({...itor(0n,128),e:-1})).toEqual({s:1,e:0n,p:64,m:1n << 63n});
  expect(rtodbl(expr(itor(0n,64)))).toBe(1);
  expect(rtodbl({s:1,e:-(1n << 60n),p:64,m:1n << 63n})).toBe(0);
  expect(() => rtodbl({s:1,e:1n << 60n,p:64,m:1n << 63n})).toThrow(PariError);
});

test('native binary splitting includes all small bases and reports rounded-ratio overflow', () => {
  const A = abpq_init(4);
  for (let i=0;i<5;i++) {A.a[i]=1n;A.b[i]=1n;A.p[i]=1n;A.q[i]=2n;}
  for (const n of [1,2,3,4]) {
    const R=abpq_sum(0,n,A);
    expect(R).toEqual({P:1n,Q:1n << BigInt(n),B:1n,T:(1n << BigInt(n))-1n});
  }
  expect(() => atanhuu((1n << 63n)-1n,(1n << 63n)+1n,64)).toThrow(new PariError('overflow in atanhuu'));
  expect(mplog2(64)).toEqual({s:1,e:-1,p:64,m:12786308645202655660n});
});
