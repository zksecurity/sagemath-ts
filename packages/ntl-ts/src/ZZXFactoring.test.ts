import {expect,test} from 'bun:test';
import {ZZX_GCD,ZZX_SquareFreeDecomp} from './index.js';
test('native integer GCD retains content and positive leading coefficient',()=>{
  expect(ZZX_GCD([-2n,0n,2n],[2n,-4n,2n])).toEqual([-2n,2n]);
  expect(ZZX_GCD([],[-3n,0n,-9n])).toEqual([3n,0n,9n]);
});
test('NTL primitive squarefree decomposition retains nonmonic integer factors',()=>{
  expect(ZZX_SquareFreeDecomp([0n,1n,2n,1n])).toEqual([[[0n,1n],1],[[1n,1n],2]]);
  expect(ZZX_SquareFreeDecomp([1n,4n,4n])).toEqual([[[1n,2n],2]]);
});
