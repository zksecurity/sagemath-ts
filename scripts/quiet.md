# Bounded command output

Use `quiet` for noninteractive commands that might print a lot. It streams stdout
and stderr into a private OS temporary log and prints at most **6,000 characters**
by default, including its status and log path. Use Bun’s `--silent` launcher flag
(as below) so Bun itself does not echo a potentially large argument list. A single giant line is clipped too.
The limit is characters, not a tokenizer-specific token count or a billing estimate.

```sh
bun --silent run quiet -- bun test tests/property/seeded.test.ts
bun --silent run quiet -- bun run test:property --case arith --runs 10
bun --silent run quiet -- bun run typecheck
bun --silent run quiet --max-chars 2500 --timeout 120 -- bun test tests/quiet.test.ts
bun --silent run quiet -- git show --stat HEAD
```

The wrapper passes arguments directly to the executable without a shell. For a
pipeline, explicitly invoke a shell with correctly quoted arguments. Captured
output is not a terminal; interactive prompts and programs requiring a TTY should
be run directly. Existing commands and test assertions are unchanged.

The summary prioritizes replay settings, reported test totals, the first few
diagnostics with surrounding records, and the end of the output. It strips ANSI
formatting for display. Detection is heuristic: the **exit status** determines
whether the command succeeded, and the full log remains the source of truth.
Failure-marker counts are diagnostic-line counts, not counts of failing tests.

The command's exit code is preserved. A missing executable returns 127; other
spawn failures return 126; a wrapper/log-write failure is nonzero. `--timeout`
returns 124. Interrupts are forwarded; on macOS/Linux, cancellation targets the
process group and escalates to SIGKILL after two seconds if necessary. On Windows,
only the immediate child is terminated.

Read more without dumping the whole log:

```sh
bun --silent run quiet show /path/printed/by/the/run/output.log --match 'error TS'
bun --silent run quiet show /path/printed/by/the/run/output.log --from 120 --lines 12
bun --silent run quiet show /path/to/existing/large.json --max-chars 1500
```

`show` also has a hard character limit. Matching is literal and case-sensitive;
it searches the entire record, even inside a giant clipped line. `--from` and the
`L` labels count newline records and carriage-return progress updates (CRLF counts
once). Increase `--from` to page through matches. `--lines` is capped at 100;
`--max-chars` accepts 1,000–20,000. Inspection streams the file with bounded memory.

Logs contain raw output, including terminal escapes, in the order received from
the two pipes. They live under `sagemath-quiet-*` in the OS temporary directory,
never in the repository. Logs are not size-capped or automatically deleted: remove
the printed run directory once its diagnostics are no longer needed. Avoid
copying passing output into tracked files. Follow printed property replay paths
to retain a failing seed or minimized input, rather than the output transcript.
