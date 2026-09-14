/** Error surfaced by the PARI Gen boundary, matching cypari2's exception name. */
export class PariError extends Error {
  constructor(message: string) {
    super(message);
    this.name = 'PariError';
  }
}
