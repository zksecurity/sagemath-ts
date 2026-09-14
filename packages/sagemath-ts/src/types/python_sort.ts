/** CPython 3.12.5 Objects/listobject.c: stable adaptive powersort.
 * Source: https://github.com/python/cpython/blob/v3.12.5/Objects/listobject.c
 * Adapted under the Python Software Foundation license; see CPYTHON-LICENSE.txt.
 * @see Deviation: Complex display ordering and Python sorting
 * The comparison sequence is observable for unordered values such as NaN.
 */
export function sorted<T>(input: readonly T[], lt: (a: T, b: T) => boolean): T[] {
  const values = input.slice(),
    size = values.length;
  let minGallop = 7;
  const stack: { start: number; length: number; power: number }[] = [];
  function gallopLeft(key: T, a: T[], start: number, n: number, hint: number): number {
    let last = 0,
      offset = 1;
    if (lt(a[start + hint]!, key)) {
      const max = n - hint;
      while (offset < max && lt(a[start + hint + offset]!, key)) {
        last = offset;
        offset = 2 * offset + 1;
      }
      offset = Math.min(offset, max);
      last += hint;
      offset += hint;
    } else {
      const max = hint + 1;
      while (offset < max && !lt(a[start + hint - offset]!, key)) {
        last = offset;
        offset = 2 * offset + 1;
      }
      offset = Math.min(offset, max);
      [last, offset] = [hint - offset, hint - last];
    }
    last++;
    while (last < offset) {
      const middle = last + Math.floor((offset - last) / 2);
      if (lt(a[start + middle]!, key)) last = middle + 1;
      else offset = middle;
    }
    return offset;
  }
  function gallopRight(key: T, a: T[], start: number, n: number, hint: number): number {
    let last = 0,
      offset = 1;
    if (lt(key, a[start + hint]!)) {
      const max = hint + 1;
      while (offset < max && lt(key, a[start + hint - offset]!)) {
        last = offset;
        offset = 2 * offset + 1;
      }
      offset = Math.min(offset, max);
      [last, offset] = [hint - offset, hint - last];
    } else {
      const max = n - hint;
      while (offset < max && !lt(key, a[start + hint + offset]!)) {
        last = offset;
        offset = 2 * offset + 1;
      }
      offset = Math.min(offset, max);
      last += hint;
      offset += hint;
    }
    last++;
    while (last < offset) {
      const middle = last + Math.floor((offset - last) / 2);
      if (lt(key, a[start + middle]!)) offset = middle;
      else last = middle + 1;
    }
    return offset;
  }
  function mergeLo(a: number, na: number, b: number, nb: number): void {
    const temp = values.slice(a, a + na);
    let i = 0,
      j = b,
      dest = a,
      mg = minGallop;
    const finishA = () => {
      for (let k = 0; k < na; k++) values[dest + k] = temp[i + k]!;
    };
    const finishB = () => {
      values.copyWithin(dest, j, j + nb);
      values[dest + nb] = temp[i]!;
    };
    values[dest++] = values[j++]!;
    nb--;
    if (nb === 0) return finishA();
    if (na === 1) return finishB();
    for (;;) {
      let acount = 0,
        bcount = 0;
      for (;;) {
        if (lt(values[j]!, temp[i]!)) {
          values[dest++] = values[j++]!;
          bcount++;
          acount = 0;
          nb--;
          if (nb === 0) return finishA();
          if (bcount >= mg) break;
        } else {
          values[dest++] = temp[i++]!;
          acount++;
          bcount = 0;
          na--;
          if (na === 1) return finishB();
          if (acount >= mg) break;
        }
      }
      mg++;
      do {
        if (mg > 1) mg--;
        minGallop = mg;
        acount = gallopRight(values[j]!, temp, i, na, 0);
        for (let k = 0; k < acount; k++) values[dest + k] = temp[i + k]!;
        dest += acount;
        i += acount;
        na -= acount;
        if (na === 1) return finishB();
        if (na === 0) return;
        values[dest++] = values[j++]!;
        nb--;
        if (nb === 0) return finishA();
        bcount = gallopLeft(temp[i]!, values, j, nb, 0);
        values.copyWithin(dest, j, j + bcount);
        dest += bcount;
        j += bcount;
        nb -= bcount;
        if (nb === 0) return finishA();
        values[dest++] = temp[i++]!;
        na--;
        if (na === 1) return finishB();
      } while (acount >= 7 || bcount >= 7);
      minGallop = ++mg;
    }
  }
  function mergeHi(a: number, na: number, b: number, nb: number): void {
    const temp = values.slice(b, b + nb);
    let i = a + na - 1,
      j = nb - 1,
      dest = b + nb - 1,
      mg = minGallop;
    const finishB = () => {
      for (let k = 0; k < nb; k++) values[dest - nb + 1 + k] = temp[k]!;
    };
    const finishA = () => {
      values.copyWithin(dest - na + 1, i - na + 1, i + 1);
      values[dest - na] = temp[j]!;
    };
    values[dest--] = values[i--]!;
    na--;
    if (na === 0) return finishB();
    if (nb === 1) return finishA();
    for (;;) {
      let acount = 0,
        bcount = 0;
      for (;;) {
        if (lt(temp[j]!, values[i]!)) {
          values[dest--] = values[i--]!;
          acount++;
          bcount = 0;
          na--;
          if (na === 0) return finishB();
          if (acount >= mg) break;
        } else {
          values[dest--] = temp[j--]!;
          bcount++;
          acount = 0;
          nb--;
          if (nb === 1) return finishA();
          if (bcount >= mg) break;
        }
      }
      mg++;
      do {
        if (mg > 1) mg--;
        minGallop = mg;
        acount = na - gallopRight(temp[j]!, values, a, na, na - 1);
        values.copyWithin(dest - acount + 1, i - acount + 1, i + 1);
        dest -= acount;
        i -= acount;
        na -= acount;
        if (na === 0) return finishB();
        values[dest--] = temp[j--]!;
        nb--;
        if (nb === 1) return finishA();
        bcount = nb - gallopLeft(values[i]!, temp, 0, nb, nb - 1);
        for (let k = 0; k < bcount; k++) values[dest - bcount + 1 + k] = temp[j - bcount + 1 + k]!;
        dest -= bcount;
        j -= bcount;
        nb -= bcount;
        if (nb === 1) return finishA();
        if (nb === 0) return;
        values[dest--] = values[i--]!;
        na--;
        if (na === 0) return finishB();
      } while (acount >= 7 || bcount >= 7);
      minGallop = ++mg;
    }
  }
  function mergeAt(index: number): void {
    const left = stack[index]!,
      right = stack[index + 1]!;
    let a = left.start,
      na = left.length,
      nb = right.length;
    const b = right.start;
    left.length += nb;
    stack.splice(index + 1, 1);
    const skip = gallopRight(values[b]!, values, a, na, 0);
    a += skip;
    na -= skip;
    if (na === 0) return;
    nb = gallopLeft(values[a + na - 1]!, values, b, nb, nb - 1);
    if (nb === 0) return;
    if (na <= nb) mergeLo(a, na, b, nb);
    else mergeHi(a, na, b, nb);
  }
  function powerloop(start: number, n1: number, n2: number): number {
    let a = 2 * start + n1,
      b = a + n1 + n2,
      power = 0;
    for (;;) {
      power++;
      if (a >= size) {
        a -= size;
        b -= size;
      } else if (b >= size) return power;
      a *= 2;
      b *= 2;
    }
  }
  let minrun = size,
    extra = 0;
  while (minrun >= 64) {
    extra = extra || minrun % 2;
    minrun = Math.floor(minrun / 2);
  }
  minrun += extra;
  for (let start = 0; start < size; ) {
    let end = start + 1;
    if (end < size) {
      const descending = lt(values[end]!, values[start]!);
      end++;
      while (end < size && lt(values[end]!, values[end - 1]!) === descending) end++;
      if (descending) {
        for (let i = start, j = end - 1; i < j; i++, j--)
          [values[i], values[j]] = [values[j]!, values[i]!];
      }
    }
    if (end - start < minrun) {
      const stop = Math.min(size, start + minrun);
      for (; end < stop; end++) {
        const pivot = values[end]!;
        let left = start,
          right = end;
        while (left < right) {
          const middle = left + Math.floor((right - left) / 2);
          if (lt(pivot, values[middle]!)) right = middle;
          else left = middle + 1;
        }
        values.copyWithin(left + 1, left, end);
        values[left] = pivot;
      }
    }
    const length = end - start;
    if (stack.length) {
      const top = stack[stack.length - 1]!,
        power = powerloop(top.start, top.length, length);
      while (stack.length > 1 && stack[stack.length - 2]!.power > power) mergeAt(stack.length - 2);
      stack[stack.length - 1]!.power = power;
    }
    stack.push({ start, length, power: 0 });
    start = end;
  }
  while (stack.length > 1) {
    let i = stack.length - 2;
    if (i > 0 && stack[i - 1]!.length < stack[i + 1]!.length) i--;
    mergeAt(i);
  }
  return values;
}
