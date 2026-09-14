#!/usr/bin/env bun
/** Capture noisy commands without buffering their output or flooding the terminal. */
import { spawn } from 'node:child_process';
import { closeSync, createReadStream, mkdtempSync, openSync, writeSync } from 'node:fs';
import { constants, tmpdir } from 'node:os';
import { basename, join } from 'node:path';
import { StringDecoder } from 'node:string_decoder';
import { stripVTControlCharacters } from 'node:util';

const DEFAULT_LIMIT = 6000;
const LINE_LIMIT = 320;
const clip = (text: string, limit: number) => limit <= 0 ? '' : text.length <= limit ? text : text.slice(0, limit - 1) + '…';
const plain = (text: string) => stripVTControlCharacters(text).replace(/[\x00-\x08\x0b-\x1f\x7f]/g, '');
type Line = { number: number; text: string; matched: boolean };

/** Retain only a line's prefix/suffix, even when a process never emits a newline. */
class Lines {
  private decoder = new StringDecoder('utf8');
  private prefix = '';
  private suffix = '';
  private overlap = '';
  private length = 0;
  private matched = false;
  private matchExcerpt = '';
  private previousCR = false;
  count = 0;
  constructor(private emit: (line: Line) => void, private match?: string) {}

  push(bytes: Buffer) { this.consume(this.decoder.write(bytes)); }
  end() { this.consume(this.decoder.end()); if (this.length) this.finish(); }

  private consume(text: string) {
    let start = 0;
    for (let i = 0; i < text.length; i++) {
      if (text[i] !== '\n' && text[i] !== '\r') continue;
      this.append(text.slice(start, i));
      // Treat carriage-return progress updates as records too; ignore the LF of CRLF.
      if (!(text[i] === '\n' && this.previousCR)) this.finish();
      this.previousCR = text[i] === '\r';
      start = i + 1;
    }
    this.append(text.slice(start));
  }

  private append(text: string) {
    if (text.length) this.previousCR = false;
    this.length += text.length;
    if (this.prefix.length < LINE_LIMIT) this.prefix += text.slice(0, LINE_LIMIT - this.prefix.length);
    this.suffix = (this.suffix + text.slice(-80)).slice(-80);
    if (this.match !== undefined) {
      const search = this.overlap + text;
      const index = search.indexOf(this.match);
      if (!this.matched && index >= 0) {
        this.matched = true;
        this.matchExcerpt = search.slice(Math.max(0, index - 30), index + this.match.length + 30);
      }
      this.overlap = this.match.length > 1 ? search.slice(-(this.match.length - 1)) : '';
    }
  }

  private finish() {
    const text = this.length > LINE_LIMIT
      ? (this.matchExcerpt ? `… match: ${this.matchExcerpt}` : this.prefix.slice(0, LINE_LIMIT - 110)) + ` … [${this.length} chars] … ` + this.suffix
      : this.prefix;
    this.emit({ number: ++this.count, text: clip(plain(text), LINE_LIMIT), matched: this.matched });
    this.prefix = this.suffix = this.overlap = '';
    this.matchExcerpt = '';
    this.length = 0; this.matched = false;
  }
}

class Summary {
  readonly head: Line[] = [];
  readonly tail: Line[] = [];
  readonly replay: Line[] = [];
  readonly totals: Line[] = [];
  readonly failures: Line[] = [];
  private previous: Line[] = [];
  private context = 0;
  failureMarkers = 0;

  add(line: Line) {
    if (!line.text.trim()) return;
    if (this.head.length < 3) this.head.push(line);
    this.tail.push(line); if (this.tail.length > 8) this.tail.shift();
    if (/\b(?:Property seed:|Native property seed:|Replay failing|Replay failed|SAGEMATH_TEST_SEED=|seed:.*path:)/i.test(line.text)) {
      if (this.replay.length < 6) this.replay.push(line);
    }
    if (/^\s*(?:\d+ (?:pass|fail|skip|todo)\b|Ran \d+ tests|Total:)|Exited with code \d+|Test (?:Files|Suites).*\d/.test(line.text)) {
      this.totals.push(line); if (this.totals.length > 8) this.totals.shift();
    }
    const failure = /\(fail\)|^\s*FAIL\b|\berror TS\d+\b|\b(?:Error|TypeError|RangeError|AssertionError):|^\s*(?:error:|fatal:)|timed out/i.test(line.text);
    if (failure) {
      this.failureMarkers++;
      if (this.failureMarkers <= 3) {
        for (const previous of this.previous) this.remember(previous);
        this.context = 5;
      }
    }
    if (this.context > 0) { this.remember(line); this.context--; }
    this.previous.push(line); if (this.previous.length > 2) this.previous.shift();
  }

  private remember(line: Line) {
    if (!this.failures.some(saved => saved.number === line.number)) this.failures.push(line);
  }

  render(limit: number): string {
    const sections: string[] = [];
    const seen = new Set<number>();
    const add = (title: string, rows: Line[], budget: number) => {
      const unique = rows.filter(row => !seen.has(row.number));
      if (!unique.length) return;
      sections.push(title + '\n' + clip(unique.map(row => {
        seen.add(row.number); return `L${row.number}: ${row.text}`;
      }).join('\n'), budget));
    };
    add('Replay:', this.replay, 1000);
    add('Reported totals:', this.totals, 800);
    add(`Diagnostics (first excerpts; ${this.failureMarkers} matching lines, not a test count):`, this.failures, 2100);
    if (!this.totals.length && !this.failures.length) add('Start:', this.head, 600);
    add('Tail:', this.tail, 1200);
    return clip(sections.join('\n\n'), limit);
  }
}

function positive(value: string | undefined, name: string, max: number): number {
  const number = Number(value);
  if (!Number.isSafeInteger(number) || number < 1 || number > max) throw Error(`${name} must be an integer in 1..${max}`);
  return number;
}

async function run(command: string[], limit: number, timeout?: number): Promise<number> {
  if (!command.length) throw Error('Expected a command after --');
  const directory = mkdtempSync(join(tmpdir(), 'sagemath-quiet-'));
  const log = join(directory, 'output.log');
  const fd = openSync(log, 'wx', 0o600);
  const start = Date.now();
  const intro = clip(`Running ${basename(command[0]!)} (${command.length - 1} arguments). Log: ${log}\n`, 500);
  process.stdout.write(intro);
  const summary = new Summary();
  const lines = new Lines(line => summary.add(line));
  let bytes = 0;
  let interrupted: NodeJS.Signals | undefined;
  let timedOut = false;
  let logError: string | undefined;
  const child = spawn(command[0]!, command.slice(1), {
    shell: false, stdio: ['inherit', 'pipe', 'pipe'], detached: process.platform !== 'win32',
  });
  const kill = (signal: NodeJS.Signals) => {
    try {
      if (process.platform !== 'win32' && child.pid) process.kill(-child.pid, signal);
      else child.kill(signal);
    } catch { /* The child may have exited between the signal and delivery. */ }
  };
  let escalation: ReturnType<typeof setTimeout> | undefined;
  const stop = (signal: NodeJS.Signals) => {
    kill(signal);
    escalation ??= setTimeout(() => kill('SIGKILL'), 2000);
  };
  const onInterrupt = () => { interrupted = 'SIGINT'; stop('SIGINT'); };
  const onTerminate = () => { interrupted = 'SIGTERM'; stop('SIGTERM'); };
  process.on('SIGINT', onInterrupt); process.on('SIGTERM', onTerminate);
  const timer = timeout === undefined ? undefined : setTimeout(() => { timedOut = true; stop('SIGTERM'); }, timeout * 1000);
  const capture = (chunk: Buffer) => {
    bytes += chunk.length;
    if (!logError) try {
      let offset = 0;
      while (offset < chunk.length) offset += writeSync(fd, chunk, offset, chunk.length - offset);
    } catch (error) { logError = String(error); stop('SIGTERM'); }
    lines.push(chunk);
  };
  child.stdout!.on('data', capture); child.stderr!.on('data', capture);
  let spawnError: NodeJS.ErrnoException | undefined;
  child.on('error', error => { spawnError = error; });
  const [code, signal] = await new Promise<[number | null, NodeJS.Signals | null]>(resolve => {
    child.on('close', (code, signal) => resolve([code, signal]));
  });
  // Descendants with redirected stdio may outlive the parent without holding its
  // pipes open. Do not cancel escalation and leave them running after cancellation.
  if (timedOut || interrupted || logError) kill('SIGKILL');
  clearTimeout(timer); clearTimeout(escalation);
  process.off('SIGINT', onInterrupt); process.off('SIGTERM', onTerminate);
  lines.end(); closeSync(fd);
  const exit = timedOut ? 124 : logError ? 1 : spawnError ? (spawnError.code === 'ENOENT' ? 127 : 126)
    : interrupted ? 128 + constants.signals[interrupted] : code ?? 128 + (signal ? constants.signals[signal] : 1);
  const footer = `\nFull log: ${log}\nInspect: bun --silent run quiet show ${JSON.stringify(log)} --match "error"\n`;
  const status = `Exit: ${exit}${timedOut ? ' (timeout)' : signal ? ` (${signal})` : ''}; ${((Date.now() - start) / 1000).toFixed(2)}s; captured ${bytes} bytes / ${lines.count} records.\n`;
  const diagnostic = spawnError || logError ? clip(plain(String(spawnError ?? logError)), 400) + '\n' : '';
  const budget = limit - intro.length - footer.length - status.length - diagnostic.length;
  process.stdout.write(clip(status + diagnostic + summary.render(Math.max(0, budget)) + footer, limit - intro.length));
  return exit;
}

async function show(path: string, options: string[]): Promise<void> {
  let limit = DEFAULT_LIMIT, from = 1, count = 20, match: string | undefined;
  for (let i = 0; i < options.length; i++) {
    const option = options[i];
    if (option === '--max-chars') limit = positive(options[++i], option, 20000);
    else if (option === '--from') from = positive(options[++i], option, Number.MAX_SAFE_INTEGER);
    else if (option === '--lines') count = positive(options[++i], option, 100);
    else if (option === '--match') {
      match = options[++i];
      if (!match || match.length > 200) throw Error('--match requires 1..200 characters (literal, case-sensitive)');
    } else throw Error(`Unknown option: ${option}`);
  }
  if (limit < 1000) throw Error('--max-chars must be at least 1000');
  const rows: Line[] = [];
  const lines = new Lines(line => {
    if (rows.length < count && line.number >= from && (match === undefined || line.matched)) rows.push(line);
  }, match);
  for await (const chunk of createReadStream(path)) lines.push(chunk as Buffer);
  lines.end();
  const footer = `\nScanned ${lines.count} records; selected ${rows.length}. Long records are clipped.\n`;
  process.stdout.write(clip(rows.map(row => `L${row.number}: ${row.text}`).join('\n') || 'No matching records.', limit - footer.length) + footer);
}

async function main(): Promise<void> {
  const args = process.argv.slice(2);
  if (!args.length || args[0] === '--help') {
    console.log('Usage: bun --silent run quiet [--max-chars 6000] [--timeout SECONDS] -- COMMAND [ARGS...]\n       bun --silent run quiet show LOG [--match TEXT] [--from RECORD] [--lines 20] [--max-chars 6000]\nRuns without a shell; keeps full output in an OS temporary log; preserves exit codes.\nSummaries are heuristic, capped at 6000 characters by default. Timeout exits 124.');
    return;
  }
  if (args[0] === 'show') {
    if (!args[1]) throw Error('Expected a log path');
    await show(args[1], args.slice(2)); return;
  }
  let limit = DEFAULT_LIMIT, timeout: number | undefined, index = 0;
  for (; index < args.length; index++) {
    const option = args[index];
    if (option === '--') { index++; break; }
    if (!option!.startsWith('-')) break; // Bun may consume the first -- itself.
    if (option === '--max-chars') limit = positive(args[++index], option, 20000);
    else if (option === '--timeout') timeout = positive(args[++index], option, 2_147_483);
    else throw Error(`Unknown option: ${option}; put -- before the command`);
  }
  if (limit < 1000) throw Error('--max-chars must be at least 1000');
  process.exitCode = await run(args.slice(index), limit, timeout);
}

if (import.meta.main) main().catch(error => { console.error(clip(plain(String(error)), 500)); process.exitCode = 2; });
