/**
 * @see Deviation: Polynomial String Parsing
 * @module sage/misc/parser
 * Arithmetic expression parser from reference/sage/src/sage/misc/parser.pyx.
 * Operator callbacks replace Python's overloaded arithmetic; no eval is used.
 */
import { NotImplementedError, OverflowError, ValueError } from '../errors.js';

export const INT = 128,
  FLOAT = 129,
  NAME = 130,
  EOS = 131,
  ERROR = 132;
export const LESS_EQ = 133,
  GREATER_EQ = 134,
  NOT_EQ = 135,
  MATRIX = 136;
export function token_to_str(token: number): string {
  if (token < -2147483648 || token > 2147483647)
    throw new OverflowError('value too large to convert to int');
  if (token < 0 || token > 0x10ffff) throw new ValueError('chr() arg not in range(0x110000)');
  return (
    (
      {
        128: 'INT',
        129: 'FLOAT',
        130: 'NAME',
        131: 'EOS',
        132: 'ERROR',
        133: 'LESS_EQ',
        134: 'GREATER_EQ',
        135: 'NOT_EQ',
        136: 'MATRIX',
      } as Record<number, string>
    )[token] ?? String.fromCodePoint(token)
  );
}
// Python's Unicode whitespace set; unlike JS whitespace it excludes BOM.
const whitespace =
  // biome-ignore lint/suspicious/noControlCharactersInRegex: Python explicitly classifies these controls as whitespace.
  /[\u0009-\u000d\u001c-\u0020\u0085\u00a0\u1680\u2000-\u200a\u2028\u2029\u202f\u205f\u3000]/u;
const alphanumeric = /[\p{L}\p{N}_]/u;
// Python isdigit also includes the Unicode Digit characters outside Nd.
// Enumerated with Sage's unicodedata; Numeric-only characters remain names.
const digit =
  /[\p{Nd}²³¹፩፪፫፬፭፮፯፰፱᧚⁰⁴⁵⁶⁷⁸⁹₀₁₂₃₄₅₆₇₈₉①②③④⑤⑥⑦⑧⑨⑴⑵⑶⑷⑸⑹⑺⑻⑼⒈⒉⒊⒋⒌⒍⒎⒏⒐⓪⓵⓶⓷⓸⓹⓺⓻⓼⓽⓿❶❷❸❹❺❻❼❽❾➀➁➂➃➄➅➆➇➈➊➋➌➍➎➏➐➑➒𐩀𐩁𐩂𐩃𐹠𐹡𐹢𐹣𐹤𐹥𐹦𐹧𐹨𑁒𑁓𑁔𑁕𑁖𑁗𑁘𑁙𑁚🄀🄁🄂🄃🄄🄅🄆🄇🄈🄉🄊]/u;
export class Tokenizer {
  readonly s: string;
  readonly chars: string[];
  pos = 0;
  last_pos = 0;
  token = 0;
  constructor(s: string) {
    this.s = s;
    this.chars = Array.from(s);
  }
  test(): string[] {
    const result: string[] = [];
    for (let token = this.next(); token !== EOS; token = this.next())
      result.push(
        [INT, FLOAT, NAME].includes(token)
          ? token_to_str(token) + '(' + this.last_token_string() + ')'
          : token_to_str(token)
      );
    return result;
  }
  reset(pos = 0): void {
    this.pos = this.last_pos = pos;
  }
  private find(): number {
    while (this.pos < this.chars.length && whitespace.test(this.chars[this.pos]!)) this.pos++;
    if (this.pos === this.chars.length) return EOS;
    const start = this.pos;
    const c = this.chars[start]!;
    const next = this.chars[start + 1];
    const pairs: Record<string, number> = {
      '<=': LESS_EQ,
      '>=': GREATER_EQ,
      '!=': NOT_EQ,
      '==': 61,
      '**': 94,
    };
    const pair = pairs[c + (next ?? '')];
    if (pair !== undefined) {
      this.pos += 2;
      return pair;
    }
    if ('+-*/^()=><,[]{}!'.includes(c)) {
      this.pos++;
      return c.codePointAt(0)!;
    }
    if (digit.test(c) || c === '.') {
      let kind = INT;
      let seenExp = false;
      let seenDecimal = false;
      while (this.pos < this.chars.length) {
        const ch = this.chars[this.pos]!;
        if (digit.test(ch)) {
        } else if (ch === '.') {
          if (seenExp || seenDecimal) break;
          kind = FLOAT;
          seenDecimal = true;
        } else if (ch === 'e' || ch === 'E') {
          if (seenExp) break;
          kind = FLOAT;
          seenExp = true;
        } else if (ch === '+' || ch === '-') {
          if (!seenExp || !['e', 'E'].includes(this.chars[this.pos - 1]!)) break;
        } else break;
        this.pos++;
      }
      return kind;
    }
    if (alphanumeric.test(c)) {
      while (this.pos < this.chars.length && alphanumeric.test(this.chars[this.pos]!)) this.pos++;
      return this.chars.slice(start, this.pos).join('') === 'matrix' ? MATRIX : NAME;
    }
    this.pos++;
    return ERROR;
  }
  next(): number {
    while (this.pos < this.chars.length && whitespace.test(this.chars[this.pos]!)) this.pos++;
    this.last_pos = this.pos;
    this.token = this.find();
    return this.token;
  }
  last(): number {
    return this.token;
  }
  peek(): number {
    const pos = this.pos;
    const token = this.find();
    this.pos = pos;
    return token;
  }
  backtrack(): boolean {
    if (this.pos === this.last_pos && this.token !== EOS)
      throw new NotImplementedError('Can only backtrack once.');
    this.pos = this.last_pos;
    this.token = 0;
    return false;
  }
  last_token_string(): string {
    return this.chars.slice(this.last_pos, this.pos).join('');
  }
}
export interface ParserOperations<T> {
  binary(op: string, left: T, right: T): T;
  unary(op: string, value: T): T;
}
export class LookupNameMaker<T> {
  constructor(
    private names: Record<string, T>,
    private fallback?: (name: string) => T
  ) {}
  set_names(names: Record<string, T>): void {
    this.names = names;
  }
  __call__(name: string): T {
    if (Object.hasOwn(this.names, name)) return this.names[name]!;
    if (this.fallback) return this.fallback(name);
    const error = new Error("Unknown variable: '" + name + "'");
    error.name = 'NameError';
    throw error;
  }
}
export class Parser<T> {
  private variable_constructor: (name: string) => T;
  private callable_constructor: (name: string) => (...args: T[]) => T;
  constructor(
    private make_int: (s: string) => T,
    private make_float: (s: string) => T,
    make_var: ((s: string) => T) | Record<string, T>,
    make_function: ((s: string) => (...args: T[]) => T) | Record<string, (...args: T[]) => T> = {},
    private implicit_multiplication = true,
    private operations: ParserOperations<T>
  ) {
    const variables = typeof make_var === 'function' ? undefined : new LookupNameMaker(make_var);
    const functions =
      typeof make_function === 'function' ? undefined : new LookupNameMaker(make_function);
    this.variable_constructor =
      typeof make_var === 'function' ? make_var : (name) => variables!.__call__(name);
    this.callable_constructor =
      typeof make_function === 'function' ? make_function : (name) => functions!.__call__(name);
  }
  _variable_constructor(): (name: string) => T {
    return this.variable_constructor;
  }
  _callable_constructor(): (name: string) => (...args: T[]) => T {
    return this.callable_constructor;
  }
  parse_sequence(_source: string): unknown {
    throw new NotImplementedError('SAGE_NOT_IMPLEMENTED: Parser.parse_sequence');
  }
  p_matrix(_tokens: Tokenizer): unknown {
    throw new NotImplementedError('SAGE_NOT_IMPLEMENTED: Parser.p_matrix');
  }
  p_sequence(_tokens: Tokenizer): unknown {
    throw new NotImplementedError('SAGE_NOT_IMPLEMENTED: Parser.p_sequence');
  }
  p_list(_tokens: Tokenizer): unknown {
    throw new NotImplementedError('SAGE_NOT_IMPLEMENTED: Parser.p_list');
  }
  p_tuple(_tokens: Tokenizer): unknown {
    throw new NotImplementedError('SAGE_NOT_IMPLEMENTED: Parser.p_tuple');
  }
  p_args(_tokens: Tokenizer): unknown {
    throw new NotImplementedError('SAGE_NOT_IMPLEMENTED: Parser.p_args');
  }
  p_arg(_tokens: Tokenizer): unknown {
    throw new NotImplementedError('SAGE_NOT_IMPLEMENTED: Parser.p_arg');
  }
  parse(s: string, accept_eqn = true): T {
    const tokens = new Tokenizer(s);
    if (tokens.peek() === MATRIX) {
      tokens.next();
      return this.p_matrix(tokens) as T;
    }
    const result = accept_eqn ? this.p_eqn(tokens) : this.p_expr(tokens);
    if (tokens.next() !== EOS) this.parse_error(tokens);
    return result;
  }
  parse_expression(s: string): T {
    return this.parse(s, false);
  }
  p_eqn(tokens: Tokenizer): T {
    const left = this.p_expr(tokens);
    const op = tokens.next();
    if ([61, NOT_EQ, 60, LESS_EQ, 62, GREATER_EQ].includes(op))
      return this.operations.binary(token_to_str(op), left, this.p_expr(tokens));
    tokens.backtrack();
    return left;
  }
  p_expr(tokens: Tokenizer): T {
    let result = this.p_term(tokens);
    let op = tokens.next();
    while (op === 43 || op === 45) {
      result = this.operations.binary(token_to_str(op), result, this.p_term(tokens));
      op = tokens.next();
    }
    tokens.backtrack();
    return result;
  }
  p_term(tokens: Tokenizer): T {
    let result = this.p_factor(tokens);
    let op = tokens.next();
    if (op === NAME && this.implicit_multiplication) {
      op = 42;
      tokens.backtrack();
    }
    while (op === 42 || op === 47) {
      result = this.operations.binary(token_to_str(op), result, this.p_factor(tokens));
      op = tokens.next();
      if (op === NAME && this.implicit_multiplication) {
        op = 42;
        tokens.backtrack();
      }
    }
    tokens.backtrack();
    return result;
  }
  p_factor(tokens: Tokenizer): T {
    const op = tokens.next();
    if (op === 43) return this.p_factor(tokens);
    if (op === 45) return this.operations.unary('-', this.p_factor(tokens));
    tokens.backtrack();
    return this.p_power(tokens);
  }
  p_power(tokens: Tokenizer): T {
    let result = this.p_atom(tokens);
    const op = tokens.next();
    if (op === 94) return this.operations.binary('^', result, this.p_factor(tokens));
    if (op === 33) {
      result = this.operations.unary('!', result);
      if (tokens.peek() === 94) {
        tokens.next();
        return this.operations.binary('^', result, this.p_factor(tokens));
      }
      return result;
    }
    tokens.backtrack();
    return result;
  }
  p_atom(tokens: Tokenizer): T {
    const token = tokens.next();
    if (token === INT) return this.make_int(tokens.last_token_string());
    if (token === FLOAT) return this.make_float(tokens.last_token_string());
    if (token === NAME) {
      const name = tokens.last_token_string();
      if (tokens.next() === 40) {
        this.callable_constructor(name);
        throw new NotImplementedError('SAGE_NOT_IMPLEMENTED: Parser function calls');
      }
      tokens.backtrack();
      return this.variable_constructor(name);
    }
    if (token === 40) {
      const result = this.p_expr(tokens);
      if (tokens.next() !== 41) this.parse_error(tokens, 'Mismatched parentheses');
      return result;
    }
    return this.parse_error(tokens);
  }
  private parse_error(tokens: Tokenizer, message = 'Malformed expression'): never {
    const error = new SyntaxError(message);
    Object.assign(error, { source: tokens.s, position: tokens.pos });
    throw error;
  }
}
