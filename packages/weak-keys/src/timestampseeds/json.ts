/** Preserve arbitrary Python JSON integers, including seeds and checkpoint indices. */
export function parse_json(text: string): unknown {
  return JSON.parse(text, (_key: string, value: unknown, context?: { source?: string }) => {
    if (
      typeof value === 'number' &&
      !Number.isSafeInteger(value) &&
      context?.source &&
      /^-?\d+$/.test(context.source)
    )
      return BigInt(context.source);
    return value;
  });
}
/** Python sort_keys=True, separators=(',', ':') for public JSONL byte compatibility. */
export function json(value: unknown): string {
  if (typeof value === 'bigint') return value.toString();
  if (Array.isArray(value)) return `[${value.map(json).join(',')}]`;
  if (value !== null && typeof value === 'object')
    return `{${Object.entries(value)
      .sort(([a], [b]) => (a < b ? -1 : a > b ? 1 : 0))
      .map(([k, v]) => `${JSON.stringify(k)}:${json(v)}`)
      .join(',')}}`;
  const result = JSON.stringify(value);
  if (result === undefined) throw new Error('Invalid JSON value');
  return result;
}
