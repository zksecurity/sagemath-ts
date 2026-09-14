/** gdisplay in PARI init.c, with the native initialized breakloop setting. */
export function pariErrorPayload(value: string, type = 't_POL'): string {
  return value.length < 1600
    ? value
    : '\n  ***  (...) Huge ' + type + ' omitted; you can access it via dbg_err()';
}
