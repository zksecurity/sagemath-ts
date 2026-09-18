Historical partial work saved when the audit paused after e34ae87.

On 2026-09-18, work.patch was applied and its is_x_coord conversion repair was
completed, including specialized GF2 constructors and focused live comparative
regressions. All 1,750 historical predicate inputs now match the live original.
Do not apply this patch again.

The four gzip files preserve the original research inputs and old/partial/native
results for the still-open lift_x promotion work. Decode only into temporary
storage when needed; do not add new bulk records. See TODO.md for current scope.
The archived partial results are historical evidence, not passing snapshots.
