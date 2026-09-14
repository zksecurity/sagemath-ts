# gsl-ts

Portable TypeScript subset of the GNU Scientific Library, used by sagemath-ts for
binary64 real integer powers. Exports `gsl_pow_int`, `gsl_sf_log`, and `gsl_sf_exp`.
This package does not implement the rest of GSL or its error-handler/status APIs.

The implementation follows GSL 2.8 `sys/pow_int.c`, `specfunc/log.c` and
`specfunc/exp.c`, from https://ftp.gnu.org/gnu/gsl/gsl-2.8.tar.gz
(SHA-256 `6a99eeed15632c6354895b1dd542ed5a855c0f15d9ad1326c6fe2b2c9e423190`).
GSL is copyright the respective GNU Scientific Library authors; this port is
licensed under GPL-3.0-or-later (see COPYING).

`gsl_pow_int(x, n)` takes a JavaScript number and a signed 32-bit integer exponent.
It preserves GSL's inversion-before-multiplication and binary64 operation order.
The log/exp functions return values with GSL's domain and normal-range guards,
using JavaScript `Math.log` and `Math.exp` for the system transcendental calls.
Their last bits can depend on the JavaScript engine and native math library.

Comparative tests execute the installed original GSL functions through SageMath:
`tests/property/cases/polynomials.cases.json` (`polynomial_real_integer_power`).
See root DEVIATIONS.md, “Polynomial Modular Powers,” for the supported domain.
