/** FLINT fmpz/get.c: conversion to binary64, truncating towards zero. */
export function fmpz_get_d(value: bigint): number {
  const magnitude = value < 0n ? -value : value;
  if (magnitude <= 9007199254740992n) return Number(value);
  const bits = magnitude.toString(2).length;
  if (bits > 1024) return value < 0n ? -Infinity : Infinity;
  const shift = bits - 53;
  const result = Number(magnitude >> BigInt(shift)) * 2 ** shift;
  return value < 0n ? -result : result;
}
