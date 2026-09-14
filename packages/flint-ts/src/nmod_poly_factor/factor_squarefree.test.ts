import {expect,test} from 'bun:test';
import {nmod_poly_factor_squarefree} from '../index.js';
test('native squarefree decomposition separates characteristic-power factors',()=>{
  expect(nmod_poly_factor_squarefree([0n,1n,0n,0n,0n,0n,1n],5n)).toEqual([[[0n,1n],1],[[1n,1n],5]]);
  expect(nmod_poly_factor_squarefree([1n,0n,0n,0n,1n],2n)).toEqual([[[1n,1n],4]]);
  expect(nmod_poly_factor_squarefree([],5n)).toEqual([]);
});
