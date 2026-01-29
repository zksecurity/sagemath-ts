let cachedSage = null;

function fallbackSage() {
  function gcd(a, b) {
    a = a < 0n ? -a : a;
    b = b < 0n ? -b : b;
    while (b !== 0n) {
      [a, b] = [b, a % b];
    }
    return a;
  }

  function xgcd(a, b) {
    let x0 = 1n;
    let y0 = 0n;
    let x1 = 0n;
    let y1 = 1n;
    while (b !== 0n) {
      const q = a / b;
      [a, b] = [b, a % b];
      [x0, x1] = [x1, x0 - q * x1];
      [y0, y1] = [y1, y0 - q * y1];
    }
    return [a, x0, y0];
  }

  function inverse_mod(a, n) {
    const [g, x] = xgcd(a, n);
    if (g !== 1n && g !== -1n) throw new Error('inverse does not exist');
    return ((x % n) + n) % n;
  }

  function power_mod(base, exp, mod) {
    let result = 1n;
    base %= mod;
    while (exp > 0n) {
      if (exp & 1n) result = (result * base) % mod;
      base = (base * base) % mod;
      exp >>= 1n;
    }
    return result;
  }

  return { gcd, inverse_mod, power_mod };
}

export function resolveSage(runtimeEl) {
  if (!cachedSage) {
    cachedSage = globalThis.Sage ?? fallbackSage();
  }
  if (runtimeEl) {
    runtimeEl.textContent = cachedSage === globalThis.Sage ? 'Runtime: sagemath-ts bundle' : 'Runtime: fallback mini-kernel';
  }
  return cachedSage;
}

function formatValue(value) {
  if (typeof value === 'bigint') return value.toString();
  if (typeof value === 'object') return JSON.stringify(value, null, 2);
  return String(value);
}

export function createRunner({ outputEl, statusEl, runtimeEl }) {
  function print(...args) {
    const line = args.map(formatValue).join(' ');
    outputEl.textContent += line + '\n';
  }

  return {
    run(code) {
      outputEl.textContent = '';
      if (statusEl) statusEl.textContent = 'Running...';

      const Sage = resolveSage(runtimeEl);
      const consoleProxy = {
        log: (...args) => print(...args),
        warn: (...args) => print('warn:', ...args),
        error: (...args) => print('error:', ...args),
      };

      try {
        const fn = new Function('Sage', 'print', 'console', code);
        fn(Sage, print, consoleProxy);
        if (statusEl) statusEl.textContent = 'Done';
      } catch (err) {
        if (statusEl) statusEl.textContent = 'Error';
        print(err instanceof Error ? err.stack || err.message : String(err));
      }
    },
  };
}
