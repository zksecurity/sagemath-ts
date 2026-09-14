import { closeSync, openSync, writeFileSync, renameSync, rmSync } from 'node:fs';
import { writeGzipChunks } from './storage.js';

/** Write transcripts without building one huge UTF-16 string in JavaScriptCore. */
export function writeJSONArray(destination: string | number, values: readonly unknown[]): void {
  const fd = typeof destination === 'string' ? openSync(destination, 'w') : destination;
  try {
    for (const chunk of jsonChunks(values)) writeFileSync(fd, chunk);
  } finally {
    if (typeof destination === 'string') closeSync(fd);
  }
}

function* jsonChunks(values: readonly unknown[]): Generator<string> {
  yield '[\n';
  for (let start = 0; start < values.length; start += 64) {
    const chunk = values
      .slice(start, start + 64)
      .map((value, offset) => {
        const json = JSON.stringify(value).replace(
          /[\u007f-\uffff]/g,
          (character) => `\\u${character.charCodeAt(0).toString(16).padStart(4, '0')}`
        );
        return (start + offset ? ',\n' : '') + json;
      })
      .join('');
    yield chunk;
  }
  yield '\n]\n';
}

export async function writeJSONArrayGzip(path: string, values: readonly unknown[]): Promise<void> {
  const temporary = `${path}.tmp`;
  try {
    await writeGzipChunks(temporary, jsonChunks(values));
    renameSync(temporary, path);
  } finally {
    rmSync(temporary, { force: true });
  }
}
