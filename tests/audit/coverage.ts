#!/usr/bin/env bun
import { isCaseFile, readText } from '../property/storage.js';
/**
 * Inventory production callables, explicit stubs, comparative areas and LCOV data.
 * Usage: bun tests/audit/coverage.ts [--lcov path ...] [--output path]
 * LCOV inputs must describe the current source; do not merge stale edited files.
 * Execution coverage is reported separately from comparative/source-audit coverage.
 */
import { execFileSync } from 'node:child_process';
import { readFileSync, writeFileSync } from 'node:fs';
import ts from 'typescript';

const args = process.argv.slice(2);
const lcov: string[] = [];
let output = 'AUDIT-COVERAGE.md';
for (let i = 0; i < args.length; i++) {
  if (args[i] === '--lcov') lcov.push(args[++i]!);
  else if (args[i] === '--output') output = args[++i]!;
  else throw new Error(`Unknown argument: ${args[i]}`);
}
const files = execFileSync('rg', ['--files', 'packages'], { encoding: 'utf8' })
  .trim()
  .split('\n')
  .filter(
    (f) =>
      f.includes('/src/') && f.endsWith('.ts') && !f.endsWith('.test.ts') && !f.endsWith('.d.ts')
  )
  .sort();
const hits = new Map<string, Map<number, number>>();
for (const path of lcov) {
  let lines: Map<number, number> | undefined;
  for (const line of readFileSync(path, 'utf8').split('\n')) {
    if (line.startsWith('SF:')) {
      const file = line.slice(3).replace(`${process.cwd()}/`, '');
      lines = hits.get(file) ?? new Map();
      hits.set(file, lines);
    } else if (line.startsWith('DA:') && lines) {
      const [n, count] = line.slice(3).split(',').map(Number);
      lines.set(n!, Math.max(lines.get(n!) ?? 0, count!));
    }
  }
}
type Callable = {
  name: string;
  line: number;
  isMember: boolean;
  explicitStub: boolean;
  partialStub: boolean;
};
const inventory = files.map((file) => {
  const source = readFileSync(file, 'utf8');
  const ast = ts.createSourceFile(file, source, ts.ScriptTarget.Latest, true);
  const callables: Callable[] = [];
  function visit(node: ts.Node, owner = '') {
    const nextOwner = ts.isClassDeclaration(node) && node.name ? node.name.text : owner;
    const publicMember =
      !ts.canHaveModifiers(node) ||
      !ts
        .getModifiers(node)
        ?.some(
          (m) =>
            m.kind === ts.SyntaxKind.PrivateKeyword || m.kind === ts.SyntaxKind.ProtectedKeyword
        );
    if (
      publicMember &&
      (ts.isFunctionDeclaration(node) ||
        ts.isMethodDeclaration(node) ||
        ts.isConstructorDeclaration(node) ||
        ts.isGetAccessorDeclaration(node) ||
        ts.isSetAccessorDeclaration(node)) &&
      node.body
    ) {
      const name = ts.isConstructorDeclaration(node)
        ? 'constructor'
        : (node.name?.getText(ast) ?? '<anonymous>');
      const isUnimplemented = (statement: ts.ThrowStatement) =>
        /NotImplementedError|(?:SAGE|PARI|FLINT|NTL)_NOT_IMPLEMENTED/.test(
          statement.expression.getText(ast)
        );
      let hasUnimplemented = false;
      function inspectThrows(child: ts.Node) {
        if (ts.isThrowStatement(child) && isUnimplemented(child)) hasUnimplemented = true;
        ts.forEachChild(child, inspectThrows);
      }
      inspectThrows(node.body);
      const explicitStub =
        node.body.statements.length === 1 &&
        ts.isThrowStatement(node.body.statements[0]!) &&
        isUnimplemented(node.body.statements[0]!);
      callables.push({
        isMember: !ts.isFunctionDeclaration(node),
        name: owner ? `${owner}.${name}` : name,
        line: ast.getLineAndCharacterOfPosition(node.getStart(ast)).line + 1,
        explicitStub,
        partialStub: !explicitStub && hasUnimplemented,
      });
    }
    ts.forEachChild(node, (child) => visit(child, nextOwner));
  }
  visit(ast);
  const lines = [...(hits.get(file)?.values() ?? [])];
  return { file, callables, lines: lines.length, hit: lines.filter((n) => n > 0).length };
});
const caseFiles = execFileSync('rg', ['--files', 'tests/property/cases'], { encoding: 'utf8' })
  .trim()
  .split('\n')
  .filter(isCaseFile)
  .sort();
const areas = caseFiles.map((file) => {
  const data = JSON.parse(readText(file)) as {
    module: string;
    cases: { function: string; seeds?: number[]; rows?: unknown[] }[];
  };
  return {
    area: data.module,
    functions: new Set(data.cases.map((c) => c.function)),
    cases: data.cases.reduce((n, c) => n + (c.rows?.length ?? c.seeds?.length ?? 1), 0),
  };
});
const integer = inventory.find((f) => f.file.endsWith('/rings/integer_ring.ts'))!;
const integerMethods = integer.callables.filter(
  (c) => c.isMember && c.name.startsWith('Integer.') && c.name !== 'Integer.constructor'
);
const integerCases = areas.find((a) => a.area === 'integers')!.functions;
const missingIntegers = integerMethods.filter(
  (c) => !integerCases.has(`integer_${c.name.slice('Integer.'.length)}`)
);
const rationalMethods = inventory
  .find((f) => f.file.endsWith('/rings/rational.ts'))!
  .callables.filter((c) => c.isMember && c.name.startsWith('Rational.'));
const rationalCases = areas.find((a) => a.area === 'rationals')!.functions;
const missingRationals = rationalMethods.filter(
  (c) => !rationalCases.has(`rational_${c.name.slice('Rational.'.length)}`)
);
const qqMethods = inventory
  .find((f) => f.file.endsWith('/rings/rational_field.ts'))!
  .callables.filter((c) => c.isMember && c.name.startsWith('RationalField.'));
const allCases = new Set(areas.flatMap((a) => [...a.functions]));
const missingQQ = qqMethods.filter(
  (c) =>
    !allCases.has(
      `qq_${c.name.slice('RationalField.'.length).replace('[Symbol.iterator]', 'iterator')}`
    )
);
const zzMethods = integer.callables.filter((c) => c.isMember && c.name.startsWith('IntegerRing.'));
const zzAdapters = new Map([
  ['__call__', 'zz_string'],
  ['random_element', 'zz_random_element_range'],
]);
const missingZZ = zzMethods.filter((c) => {
  const name = c.name.slice('IntegerRing.'.length);
  return !allCases.has(zzAdapters.get(name) ?? `zz_${name}`);
});
// Explicit adapters prevent a newly added member from looking covered merely
// because another member of its class has comparative cases.
const finiteMembers = inventory
  .filter((f) => /finite_rings\/finite_field_(extension|prime)\.ts$/.test(f.file))
  .flatMap((f) => f.callables.filter((c) => c.isMember).map((c) => ({ ...c, file: f.file })));
const finiteAdapter = (file: string, member: string): string | undefined => {
  const dot = member.indexOf('.');
  const owner = member.slice(0, dot);
  const name = member.slice(dot + 1);
  const legacy = file.endsWith('finite_field_prime.ts');
  const primeElement = owner === 'PrimeFieldElement' || (legacy && owner === 'FiniteFieldElement');
  if (primeElement) return name === 'constructor' ? 'ff_direct_element' : `ff_prime_${name}`;
  if (owner === 'PrimeField' || owner === 'FiniteFieldPrime') {
    if (name === 'constructor') return 'ff_prime_constructor';
    if (name === '__call__') return 'ff_scalar_constructor';
    if (name === 'random_element') return 'ff_random';
    if (name === 'quadratic_non_residue') return 'ff_prime_quadratic_non_residue';
    if (
      [
        'zero',
        'one',
        'gen',
        'cardinality',
        'is_field',
        'toString',
        '[Symbol.iterator]',
        'elements',
        'list',
        'multiplicative_generator',
        'primitive_element',
      ].includes(name)
    )
      return 'ff_prime_parent';
  }
  if (owner === 'FiniteFieldElement') {
    if (name === 'sqrt' || name === 'is_square') return `ff_extension_${name}`;
    return name === 'constructor' ? 'ff_ext_conversion' : `ff_ext_${name}`;
  }
  if (owner === 'FiniteFieldExtension') {
    if (name === '__call__') return 'ff_scalar_constructor';
    if (name === 'fromInteger') return 'ff_ext_fromInteger';
    if (name === 'random_element') return 'ff_random';
    if (name === '[Symbol.iterator]') return 'ff_ext_iteration';
    if (
      [
        'constructor',
        'zero',
        'one',
        'gen',
        'cardinality',
        'is_field',
        'toString',
        'elements',
        'primitiveElement',
        'multiplicative_generator',
        'primitive_element',
      ].includes(name)
    )
      return 'ff_ext_parent';
  }
  return undefined;
};
const missingFinite = finiteMembers.filter(
  (c) => !allCases.has(finiteAdapter(c.file, c.name) ?? '')
);
const parserCallables = inventory.find((f) => f.file.endsWith('/misc/parser.ts'))?.callables ?? [];
const parserImplemented = parserCallables.filter((c) => !c.explicitStub);
const parserAdapters = new Map<string, string>([
  ['token_to_str', 'parser_token_name'],
  ...['constructor', 'test', 'reset', 'next', 'last', 'peek', 'backtrack', 'last_token_string'].map(
    (name) => ['Tokenizer.' + name, 'parser_tokens'] as [string, string]
  ),
  ...['constructor', 'set_names', '__call__'].map(
    (name) => ['LookupNameMaker.' + name, 'parser_lookup'] as [string, string]
  ),
  ...[
    'constructor',
    '_variable_constructor',
    '_callable_constructor',
    'parse',
    'parse_expression',
    'p_eqn',
    'p_expr',
    'p_term',
    'p_factor',
    'p_power',
    'p_atom',
  ].map((name) => ['Parser.' + name, 'parser_ast'] as [string, string]),
]);
const missingParser = parserImplemented.filter(
  (c) => !allCases.has(parserAdapters.get(c.name) ?? '')
);
const modularMembers = inventory
  .filter((f) => /\/integer_mod(?:_ring)?\.ts$/.test(f.file))
  .flatMap((f) => f.callables.filter((c) => c.isMember));
const modularAdapters = new Map<string, string>([
  ['IntegerMod._rational_', 'modular__rational_'],
  ['IntegerMod.constructor', 'mi_constructor'],
  ['IntegerMod.log', 'im_log'],
  ...[
    'add',
    'sub',
    'mul',
    'div',
    'eq',
    'pow',
    'neg',
    'inv',
    'isZero',
    'isOne',
    'isUnit',
    'modulus',
    'lift',
    'toBigInt',
    'toString',
    'repr',
    'multiplicative_order',
  ].map((name) => ['IntegerMod.' + name, 'mi_' + name] as [string, string]),
  ['IntegerModRing.constructor', 'mi_ring_direct'],
  ['IntegerModRing.__call__', 'mi_constructor'],
  ['IntegerModRing.random_element', 'mi_random'],
  ['IntegerModRing.field', 'mi_field'],
  ['IntegerModRing.factored_order', 'mi_factored_order'],
  ['IntegerModRing._lift_residue_field_root', 'mi_residue_root_lift'],
  ['IntegerModRing._roots_univariate_polynomial', 'mi_modular_roots_hook'],
  ['IntegerModRing.multiplicative_generator', 'imr_mult_gen'],
  ...[
    'zero',
    'one',
    'gen',
    'is_field',
    'cardinality',
    '[Symbol.iterator]',
    'list',
    'units',
    'multiplicative_group_is_cyclic',
    'unit_gens',
    'toString',
  ].map((name) => ['IntegerModRing.' + name, 'mi_ring'] as [string, string]),
]);
const missingModular = modularMembers.filter(
  (c) => !allCases.has(modularAdapters.get(c.name) ?? '')
);
// Binary API names come from direct class/export declarations, excluding file-I/O stubs.
const binaryFile = 'packages/sagemath-ts/src/matrix/matrix_mod2.ts';
const binaryAst = ts.createSourceFile(
  binaryFile,
  readFileSync(binaryFile, 'utf8'),
  ts.ScriptTarget.Latest,
  true
);
const binaryCallables = inventory.find((f) => f.file === binaryFile)!.callables;
const binaryAPI: string[] = [];
for (const statement of binaryAst.statements) {
  if (ts.isClassDeclaration(statement) && statement.name?.text === 'Matrix_mod2_dense') {
    for (const member of statement.members) {
      if (
        (ts.isMethodDeclaration(member) || ts.isConstructorDeclaration(member)) &&
        member.body &&
        !member.modifiers?.some(
          (m) =>
            m.kind === ts.SyntaxKind.PrivateKeyword || m.kind === ts.SyntaxKind.ProtectedKeyword
        )
      ) {
        binaryAPI.push(
          'Matrix_mod2_dense.' +
            (ts.isConstructorDeclaration(member) ? 'constructor' : member.name.getText(binaryAst))
        );
      }
    }
  } else if (
    ts.isFunctionDeclaration(statement) &&
    statement.name &&
    statement.body &&
    statement.modifiers?.some((m) => m.kind === ts.SyntaxKind.ExportKeyword) &&
    !binaryCallables.find((c) => c.name === statement.name!.text)?.explicitStub
  ) {
    binaryAPI.push(statement.name.text);
  }
}
const binaryAdapters = new Map<string, string>();
function binaryMembers(area: string, members: string[]) {
  for (const member of members) binaryAdapters.set('Matrix_mod2_dense.' + member, area);
}
binaryMembers('binary_matrix_cache', [
  'set_immutable',
  'is_immutable',
  'is_mutable',
  '_clear_cache',
  'doubly_lexical_ordering',
]);
binaryMembers('binary_subdivision', ['subdivide', 'subdivisions', 'get_subdivisions']);
binaryMembers('binary_constructor', ['constructor', 'list']);
binaryMembers('binary_matrix_access', ['get', 'set', 'row', 'submatrix']);
binaryMembers('binary_format', ['str']);
binaryMembers('binary_columns_state', ['columns']);
binaryMembers('binary_basic', [
  'add',
  'sub',
  'neg',
  'copy',
  'determinant',
  'transpose',
  'augment',
  'toString',
]);
binaryMembers('binary_matrix_product', [
  '_multiply_classical',
  'mul',
  '_multiply_m4rm',
  '_multiply_strassen',
]);
binaryMembers('binary_inverse', ['inverse']);
binaryMembers('binary_matrix_elimination', ['echelonize', 'echelon_form', 'pivots', 'rank']);
binaryMembers('binary_random_state', ['randomize']);
binaryMembers('binary_matrix_density', ['density']);
binaryMembers('binary_solve', ['solve_right']);
binaryMembers('binary_kernel', ['right_kernel_matrix']);
binaryMembers('binary_predicates', ['is_Gamma_free', 'eq']);
binaryMembers('binary_swap_conversion', ['swap_rows', 'swap_columns']);
binaryMembers('binary_permutation', ['permute_rows', 'permute_columns']);
for (const name of ['from_png_data', 'to_png_data']) binaryAdapters.set(name, 'binary_png_data');
for (const name of ['ple', 'pluq']) binaryAdapters.set(name, 'binary_factorization');
for (const name of ['zero_matrix_gf2', 'identity_matrix_gf2', 'random_matrix_gf2'])
  binaryAdapters.set(name, 'binary_factory');
binaryAdapters.set('matrix_gf2_from_entries', 'binary_from_entries');
const missingBinary = binaryAPI.filter((name) => !allCases.has(binaryAdapters.get(name) ?? ''));
const allLines = inventory.reduce((n, f) => n + f.lines, 0);
const allHit = inventory.reduce((n, f) => n + f.hit, 0);
const pct = (hit: number, total: number) =>
  total ? `${((100 * hit) / total).toFixed(2)}%` : 'unmeasured';
const report = [
  '# Port audit coverage inventory',
  '',
  'Generated with `bun tests/audit/coverage.ts`. This is an inventory, not a claim of full behavioral equivalence.',
  '',
  `- Modular-integer/ring member dispatch coverage: ${modularMembers.length - missingModular.length}/${modularMembers.length}; factory aliases/zero/negative orders are compared; broader coercion/category domains remain open.`,
  '- Parser/tokenizer implemented callable dispatch coverage: ' +
    (parserImplemented.length - missingParser.length) +
    '/' +
    parserImplemented.length +
    '; excludes ' +
    parserCallables.filter((c) => c.explicitStub).length +
    ' throw-only stubs. Function-call branches remain incomplete.',
  `- Production TypeScript files: ${files.length}.`,
  `- Named functions, constructors and public methods with bodies: ${inventory.reduce((n, f) => n + f.callables.length, 0)} (includes internal functions).`,
  `- Explicit throw-only stubs: ${inventory.reduce((n, f) => n + f.callables.filter((c) => c.explicitStub).length, 0)}.`,
  `- Other callables containing an explicit unimplemented branch: ${inventory.reduce((n, f) => n + f.callables.filter((c) => c.partialStub).length, 0)}.`,
  `- Comparative areas: ${areas.length}; configured seeded cases: ${areas.reduce((n, a) => n + a.cases, 0)}.`,
  `- Integer method dispatch coverage: ${integerMethods.length - missingIntegers.length}/${integerMethods.length}. Adapters and input limits are explicit in the oracle and DEVIATIONS.md.`,
  `- Rational method, accessor and constructor dispatch coverage: ${rationalMethods.length - missingRationals.length}/${rationalMethods.length}. Symbolic/complex return branches remain explicit implementation gaps.`,
  `- RationalField member dispatch coverage: ${qqMethods.length - missingQQ.length}/${qqMethods.length}; random streams are compared in rand_stats.`,
  `- IntegerRing member dispatch coverage: ${zzMethods.length - missingZZ.length}/${zzMethods.length}. Integer construction is also compared through the zz_* dispatchers with both wrapper modes.`,
  `- Finite-field parent/element member dispatch coverage: ${finiteMembers.length - missingFinite.length}/${finiteMembers.length} across both prime-field implementations and the extension backend. Backend adapters and unimplemented input branches remain explicit in DEVIATIONS.md.`,
  `- Binary matrix supported API dispatch coverage: ${binaryAPI.length - missingBinary.length}/${binaryAPI.length}; two file-I/O stubs are excluded. Numeric/container and pixel-data adapters are explicit in DEVIATIONS.md.`,
  `- Instrumented source-line execution coverage: ${allHit}/${allLines} (${pct(allHit, allLines)}). LCOV does not establish branch or behavioral coverage.`,
  '',
  'The September source review and repairs are recorded in [AUDIT-2026-09.md](AUDIT-2026-09.md).',
  'Every file below remains in the audit inventory, including files absent from the supplied execution reports.',
  'Existing implementation gaps are tracked separately in DEVIATIONS.md; stubs are not treated as implemented behavior.',
  '',
  '## Comparative areas',
  '',
  '| Area | Dispatch functions | Seeded cases |',
  '|---|---:|---:|',
  ...areas.map((a) => `| ${a.area} | ${a.functions.size} | ${a.cases} |`),
  '',
  '## Production files',
  '',
  '| File | Callables | Throw-only stubs | Partial stubs | Executed lines |',
  '|---|---:|---:|---:|---:|',
  ...inventory.map(
    (f) =>
      `| [${f.file}](${f.file}) | ${f.callables.length} | ${f.callables.filter((c) => c.explicitStub).length} | ${f.callables.filter((c) => c.partialStub).length} | ${f.lines ? `${f.hit}/${f.lines} (${pct(f.hit, f.lines)})` : 'unmeasured'} |`
  ),
  '',
  '## Integer methods missing a comparative dispatcher',
  '',
  ...(missingIntegers.length
    ? missingIntegers.map((c) => `- ${c.name}`)
    : [
        'None. Method dispatch coverage is distinct from complete input-domain or algorithm coverage.',
      ]),
  '',
];
report.push(
  '## Binary matrix API missing a comparative dispatcher',
  '',
  ...(missingBinary.length
    ? missingBinary.map((name) => '- ' + name)
    : ['None. Dispatch coverage is distinct from full input-domain coverage.']),
  ''
);
report.push(
  '## Modular-integer/ring members missing a comparative dispatcher',
  '',
  ...(missingModular.length
    ? missingModular.map((c) => '- ' + c.name)
    : ['None. Broader input-domain and category gaps remain open.']),
  ''
);
report.push(
  '## Rational members missing a comparative dispatcher',
  '',
  ...(missingRationals.length
    ? missingRationals.map((c) => `- ${c.name}`)
    : ['None. Dispatch coverage does not establish full input-domain coverage.']),
  ''
);
report.push(
  '## RationalField members missing a comparative dispatcher',
  '',
  ...(missingQQ.length ? missingQQ.map((c) => `- ${c.name}`) : ['None.']),
  ''
);
report.push(
  '## IntegerRing members missing a comparative dispatcher',
  '',
  ...(missingZZ.length ? missingZZ.map((c) => `- ${c.name}`) : ['None.']),
  ''
);
report.push(
  '## Finite-field members missing a comparative dispatcher',
  '',
  ...(missingFinite.length
    ? missingFinite.map((c) => `- ${c.file}: ${c.name}`)
    : ['None. Dispatch coverage is not full input-domain coverage.']),
  ''
);
writeFileSync(output, report.join('\n'));
console.log(
  `${output}: ${files.length} files, ${allHit}/${allLines} instrumented lines, ${integerMethods.length - missingIntegers.length}/${integerMethods.length} Integer methods`
);
