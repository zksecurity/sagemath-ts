import { readFileSync } from 'node:fs';
import { createWriteStream } from 'node:fs';
import { Readable } from 'node:stream';
import { pipeline } from 'node:stream/promises';
import { createGzip, gunzipSync } from 'node:zlib';

/** Large suites and retained transcripts use gzip; small editable suites may stay plain. */
export function readText(path: string): string {
  const data = readFileSync(path);
  return (path.endsWith('.gz') ? gunzipSync(data) : data).toString('utf8');
}

export function isCaseFile(name: string): boolean {
  return /\.cases\.json(?:\.gz)?$/.test(name);
}

export function caseArea(name: string): string {
  return name.replace(/\.cases\.json(?:\.gz)?$/, '');
}

/** Streaming compression avoids joining a large transcript into one JS string. */
export async function writeGzipChunks(path: string, chunks: Iterable<string>): Promise<void> {
  await pipeline(Readable.from(chunks), createGzip({ level: 9 }), createWriteStream(path));
}
