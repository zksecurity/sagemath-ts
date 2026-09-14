import { IndexError, NotImplementedError, RuntimeError } from '../../errors.js';
import type { RingElement } from './polynomial_element.js';

type Node =
  | { kind: 'variable' }
  | { kind: 'coefficient'; index: number }
  | { kind: 'gap' }
  | { kind: 'square'; left: number }
  | { kind: 'multiply'; left: number; right: number }
  | { kind: 'affine'; left: number; right: number; index: number };

/**
 * Sage polynomial_compiled.pyx: sparse Horner structure and binary gap filling.
 * @see Deviation: Polynomial Evaluation and Composition
 */
export class CompiledPolynomialFunction<C extends RingElement> {
  private readonly nodes: Node[] = [{ kind: 'variable' }];
  private readonly order: number[] = [];
  private readonly references: number[];
  private readonly root: number;

  constructor(
    private readonly coefficients: readonly C[],
    algorithm = 'binary'
  ) {
    if (!coefficients.length) throw new IndexError('list index out of range');
    const nodes = this.nodes;
    const add = (node: Node): number => {
      nodes.push(node);
      return nodes.length - 1;
    };
    // A max heap supplies the same ordered gaps as Sage's BinaryTree in O(log n).
    const pending = [1],
      gaps = new Map<number, number>([[1, 0]]);
    const push = (value: number) => {
      let i = pending.length;
      pending.push(value);
      while (i) {
        const parent = (i - 1) >> 1;
        if (pending[parent]! >= value) break;
        pending[i] = pending[parent]!;
        i = parent;
      }
      pending[i] = value;
    };
    const pop = (): number => {
      const result = pending[0]!,
        last = pending.pop()!;
      if (pending.length) {
        let i = 0;
        while (2 * i + 1 < pending.length) {
          let child = 2 * i + 1;
          if (child + 1 < pending.length && pending[child + 1]! > pending[child]!) child++;
          if (pending[child]! <= last) break;
          pending[i] = pending[child]!;
          i = child;
        }
        pending[i] = last;
      }
      return result;
    };
    const gap = (width: number): number => {
      const existing = gaps.get(width);
      if (existing !== undefined) return existing;
      const id = add({ kind: 'gap' });
      gaps.set(width, id);
      push(width);
      return id;
    };
    let root = add({ kind: 'coefficient', index: coefficients.length - 1 }),
      width = 0;
    for (let i = coefficients.length - 2; i > 0; i--) {
      width++;
      if (!coefficients[i]!.isZero()) {
        root = add({ kind: 'affine', left: root, right: gap(width), index: i });
        width = 0;
      }
    }
    width++;
    root = coefficients[0]!.isZero()
      ? add({ kind: 'multiply', left: root, right: gap(width) })
      : add({ kind: 'affine', left: root, right: gap(width), index: 0 });
    this.root = root;
    if (pending[0]! > 1 && algorithm !== 'binary') {
      if (algorithm === 'pippenger')
        throw new NotImplementedError(
          "Implementation of Pippenger's Algorithm is not ready for prime time."
        );
      throw new RuntimeError("Method '%s' not supported.");
    }
    while (pending.length > 1) {
      const m = pop(),
        n = pending[0]!,
        k = Math.floor(m / n),
        r = m % n,
        half = Math.floor(m / 2);
      const id = gaps.get(m)!;
      const h = m % 2 === 0 && n >= half ? gaps.get(half) : undefined;
      if (h !== undefined) nodes[id] = { kind: 'square', left: h };
      else if (r > 0) nodes[id] = { kind: 'multiply', left: gap(n * k), right: gap(r) };
      else if (k % 2 === 0) nodes[id] = { kind: 'square', left: gap((n * k) / 2) };
      else nodes[id] = { kind: 'multiply', left: gap(n * (k - 1)), right: gaps.get(n)! };
    }
    // Iterative postorder avoids the JavaScript recursion limit while preserving
    // pd_eval's left-before-right order and one evaluation per shared node.
    const seen = new Set<number>(),
      stack: [number, boolean][] = [[root, false]];
    this.references = Array<number>(nodes.length).fill(0);
    while (stack.length) {
      const [id, ready] = stack.pop()!;
      if (seen.has(id)) continue;
      const node = nodes[id]!;
      if (ready) {
        seen.add(id);
        this.order.push(id);
        continue;
      }
      stack.push([id, true]);
      if ('right' in node) {
        this.references[node.right]!++;
        stack.push([node.right, false]);
      }
      if ('left' in node) {
        this.references[node.left]!++;
        stack.push([node.left, false]);
      }
    }
  }

  private run<T>(
    variable: T,
    coefficient: (i: number) => T,
    multiply: (a: T, b: T) => T,
    square: (a: T) => T,
    affine: (a: T, b: T, i: number) => T
  ): T {
    const values = Array<T | undefined>(this.nodes.length),
      references = this.references.slice();
    for (const id of this.order) {
      const node = this.nodes[id]!;
      switch (node.kind) {
        case 'variable':
          values[id] = variable;
          break;
        case 'coefficient':
          values[id] = coefficient(node.index);
          break;
        case 'square':
          values[id] = square(values[node.left]!);
          break;
        case 'multiply':
          values[id] = multiply(values[node.left]!, values[node.right]!);
          break;
        case 'affine':
          values[id] = affine(values[node.left]!, values[node.right]!, node.index);
          break;
        case 'gap':
          throw new Error('unfilled polynomial evaluation gap');
      }
      // pd_clean releases a value after its final consumer, including shared powers.
      if ('left' in node && --references[node.left]! === 0) values[node.left] = undefined;
      if ('right' in node && --references[node.right]! === 0) values[node.right] = undefined;
    }
    return values[this.root]!;
  }

  eval(x: C): C;
  eval(
    x: unknown,
    arithmetic: { multiply(a: unknown, b: unknown): unknown; add(a: unknown, b: unknown): unknown }
  ): unknown;
  eval(
    x: unknown,
    arithmetic?: { multiply(a: unknown, b: unknown): unknown; add(a: unknown, b: unknown): unknown }
  ): unknown {
    const multiply = arithmetic?.multiply ?? ((a: unknown, b: unknown) => (a as C).mul(b as C));
    const add = arithmetic?.add ?? ((a: unknown, b: unknown) => (a as C).add(b as C));
    return this.run<unknown>(
      x,
      (i) => this.coefficients[i]!,
      multiply,
      (a) => multiply(a, a),
      (a, b, i) => add(multiply(a, b), this.coefficients[i]!)
    );
  }
  __call__(x: C): C {
    return this.eval(x);
  }
  toString(): string {
    return `CompiledPolynomialFunction(${this.run(
      'x',
      (i) => `a${i}`,
      (a, b) => `(${a}*${b})`,
      (a) => `(${a})^2`,
      (a, b, i) => `(${a}*${b}+a${i})`
    )})`;
  }
}
