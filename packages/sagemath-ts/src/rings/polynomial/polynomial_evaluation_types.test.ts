import { expect, test } from 'bun:test';
import { resolve } from 'node:path';
import ts from 'typescript';
import { Matrix_mod2_dense } from '../../matrix/matrix_mod2.js';
import { GF2 } from '../finite_rings/gf2.js';
import { PolynomialRing } from './polynomial_ring.js';

test('binary matrix evaluation and one-level argument lists retain matrix values', () => {
  const f = new PolynomialRing(GF2, 'x').__call__([1n, 1n]);
  const a = new Matrix_mod2_dense(2, 2, [
    [1, 1],
    [0, 1],
  ]);
  const direct = f.evaluate(a),
    listed = f.evaluate([a]);
  expect(direct).toBeInstanceOf(Matrix_mod2_dense);
  expect(listed).toBeInstanceOf(Matrix_mod2_dense);
  expect([0, 1].map((i) => (listed as Matrix_mod2_dense).row(i))).toEqual([
    [0, 1],
    [0, 0],
  ]);
  expect(direct).not.toBe(a);
  expect([0, 1].map((i) => a.row(i))).toEqual([
    [1, 1],
    [0, 1],
  ]);
});

test('matrix, mixed and argument-list evaluation types cannot silently narrow to scalars', () => {
  const filename = resolve(import.meta.dir, '__evaluation_types__.ts');
  const source = `
import {Matrix_mod2_dense} from '../../matrix/matrix_mod2.js';
import {Matrix_modn_dense} from '../../matrix/matrix_modn.js';
import {IntegerMatrix} from '../../matrix/matrix_integer.js';
import {Matrix} from '../../matrix/matrix_generic.js';
import {GF2,type GF2Element} from '../finite_rings/gf2.js';
import {PolynomialRing} from './polynomial_ring.js';
import type {RingElement} from './polynomial_element.js';
const f=new PolynomialRing(GF2,'x').__call__([1n,1n]);
declare const a:Matrix_mod2_dense;
declare const mixed:Matrix_mod2_dense|GF2Element;
type Result=GF2Element|number|RingElement|Matrix_mod2_dense|Matrix_modn_dense|IntegerMatrix|Matrix<RingElement>;
const direct=f.evaluate(a);direct.get(0,0);
const listed=f.evaluate([a] as const);
const combined=f.evaluate(mixed);
const listedCompatible:Result=listed;
const combinedCompatible:Result=combined;
const mixedMayBeMatrix:Extract<typeof combined,Matrix_mod2_dense>=a;
const listMayBeMatrix:Extract<typeof listed,Matrix_mod2_dense>=a;
// @ts-expect-error a matrix argument cannot have a scalar-only result
const badDirect:RingElement=f.evaluate(a);
// @ts-expect-error an argument list may contain a matrix
const badList:RingElement=f.evaluate([a] as const);
// @ts-expect-error a mixed argument may select a matrix result
const badMixed:RingElement=f.evaluate(mixed);
`;
  const options: ts.CompilerOptions = {
    strict: true,
    noEmit: true,
    target: ts.ScriptTarget.ES2022,
    module: ts.ModuleKind.NodeNext,
    moduleResolution: ts.ModuleResolutionKind.NodeNext,
    skipLibCheck: true,
  };
  const host = ts.createCompilerHost(options),
    original = host.getSourceFile.bind(host);
  host.getSourceFile = (name, language, onError, fresh) =>
    name === filename
      ? ts.createSourceFile(name, source, language, true)
      : original(name, language, onError, fresh);
  const program = ts.createProgram([filename], options, host),
    file = program.getSourceFile(filename)!;
  expect(
    program
      .getSemanticDiagnostics(file)
      .map((d) => ts.flattenDiagnosticMessageText(d.messageText, '\n'))
  ).toEqual([]);
});
