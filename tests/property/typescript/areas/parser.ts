/** Direct grammar/tokenizer comparisons with Sage's misc/parser.pyx. */
import {
  LookupNameMaker,
  Parser,
  Tokenizer,
  token_to_str,
} from '../../../../packages/sagemath-ts/src/misc/parser.js';
import { ZZ } from '../../../../packages/sagemath-ts/src/rings/integer_ring.js';

function comparison(run: () => unknown): string {
  const normalized = (x: unknown): unknown =>
    Array.isArray(x) ? x.map(normalized) : typeof x === 'boolean' || x === null ? x : String(x);
  try {
    return JSON.stringify({ value: normalized(run()) });
  } catch (error) {
    return JSON.stringify({ error: (error as Error).name, message: (error as Error).message });
  }
}
const source = (codes: bigint[]) => String.fromCodePoint(...codes.map(Number));
export const functions: Record<string, (...args: never[]) => unknown> = {
  parser_tokens: (codes: bigint[], mode: bigint, position: bigint) =>
    comparison(() => {
      const t = new Tokenizer(source(codes));
      if (mode === 0n) return [t.test(), t.test()];
      if (mode === 1n) {
        const before = t.peek();
        const first = t.next();
        const last = t.last();
        const value = t.last_token_string();
        const after = t.peek();
        return [
          token_to_str(before),
          token_to_str(first),
          token_to_str(last),
          value,
          token_to_str(after),
          t.test(),
        ];
      }
      if (mode === 2n) {
        t.test();
        t.reset(Number(position));
        return t.test();
      }
      if (mode === 3n) {
        t.next();
        t.backtrack();
        t.backtrack();
        return t.test();
      }
      const first = t.next();
      const result = t.backtrack();
      const again = t.next();
      return [token_to_str(first), result, token_to_str(again), t.last_token_string(), t.test()];
    }),
  parser_ast: (codes: bigint[], mode: bigint, implicit: bigint) =>
    comparison(() => {
      const text = source(codes);
      const p = new Parser<string>(
        (s) => 'I' + ZZ.__call__(s),
        (s) => 'F' + s,
        (s) => 'V' + s,
        {},
        implicit !== 0n,
        {
          binary: (op, a, b) => '(' + op + ' ' + a + ' ' + b + ')',
          unary: (op, a) => '(' + op + ' ' + a + ')',
        }
      );
      if (mode === 8n) return p._variable_constructor()(text);
      if (mode === 9n) return p._callable_constructor()(text);
      if (mode === 0n) return p.parse(text);
      if (mode === 1n) return p.parse_expression(text);
      const t = new Tokenizer(text);
      const method = (['p_eqn', 'p_expr', 'p_term', 'p_factor', 'p_power', 'p_atom'] as const)[
        Number(mode) - 2
      ]!;
      return [p[method](t), t.test()];
    }),
  parser_lookup: (codes: bigint[], fallback: bigint) =>
    comparison(() => {
      const name = source(codes);
      const maker = new LookupNameMaker(
        { a: 'old', constructor: 'owned' },
        fallback ? (n) => 'fallback:' + n : undefined
      );
      const first = maker.__call__(name);
      maker.set_names({ a: 'new', ['__proto__']: 'ordinary' });
      return [first, maker.__call__(name)];
    }),
  parser_token_name: (n: bigint) => comparison(() => token_to_str(Number(n))),
};

functions.parser_constructed_names = (mode: bigint) =>
  comparison(() => {
    const p = new Parser<string>(
      (s) => s,
      (s) => s,
      { a: 'A' },
      (name) => (arg) => name + ':' + arg,
      true,
      {
        binary: (op, a, b) => '(' + op + ' ' + a + ' ' + b + ')',
        unary: (op, a) => '(' + op + ' ' + a + ')',
      }
    );
    return mode === 0n ? p.parse('a+a') : p._callable_constructor()('foo')('argument');
  });
