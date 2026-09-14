/** Native 64-bit FLINT random state (flint.h.in). GMP state is unused by these kernels. */
export interface flint_rand_t {
  __randval: bigint;
  __randval2: bigint;
}
export function flint_rand_init(): flint_rand_t {
  return { __randval: 13845646450878251009n, __randval2: 13142370077570254774n };
}
export function flint_rand_set_seed(state: flint_rand_t, seed1: bigint, seed2: bigint): void {
  state.__randval = BigInt.asUintN(64, seed1);
  state.__randval2 = BigInt.asUintN(64, seed2);
}
export function flint_rand_get_seed(state: flint_rand_t): [bigint, bigint] {
  return [state.__randval, state.__randval2];
}
/** flint.h.in: flint_rand_clear is a no-op for the native word generator. */
export function flint_rand_clear(_state: flint_rand_t): void {}
