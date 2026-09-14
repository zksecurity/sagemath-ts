/** PARI nffactor.c factor recombination heuristics. */
/** Maximum subset size before switching to the native lattice recombination. */
export function cmbf_maxK(nb: number): number {
  return nb > 10 ? 3 : nb - 1;
}
