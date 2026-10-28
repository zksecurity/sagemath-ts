/** Independent small-integer oracle. No imports from the screened math implementation. */
export function classify(n: number): { prime: boolean; power: boolean } {
  if (!Number.isSafeInteger(n) || n < 0 || n > 10000) throw new Error('Oracle input out of range');
  let prime = n > 1;
  for (let divisor = 2; divisor * divisor <= n; divisor++)
    if (n % divisor === 0) {
      prime = false;
      break;
    }
  let power = false;
  for (let base = 2; base * base <= n; base++) {
    let value = base * base;
    while (value < n) value *= base;
    if (value === n) {
      power = true;
      break;
    }
  }
  return { prime, power };
}
if (import.meta.main) {
  const values: number[] = await Bun.stdin.json();
  console.log(JSON.stringify(values.map(classify)));
}
