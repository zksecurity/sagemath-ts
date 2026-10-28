/** POSIX locks keep schema-1/2 checkpoints interoperable with the original CLI. Bun only. */
import { dlopen } from 'bun:ffi';
import { randomUUID } from 'node:crypto';
import { constants, closeSync, fsyncSync, openSync, renameSync, rmSync, writeSync } from 'node:fs';
import { basename, dirname, join } from 'node:path';
import { json } from './json.js';

let library: ReturnType<typeof load> | undefined;
function load() {
  if (process.platform !== 'darwin' && process.platform !== 'linux')
    throw new Error('Timestamp checkpoint locking requires Linux or macOS');
  return dlopen(process.platform === 'darwin' ? '/usr/lib/libSystem.B.dylib' : 'libc.so.6', {
    flock: { args: ['i32', 'i32'], returns: 'i32' },
  });
}
export function lock_fd(fd: number, message = 'Output is in use by another process'): void {
  library ??= load();
  if (library.symbols.flock(fd, 2 | 4) !== 0) throw new Error(message);
}
export function checkpoint_lock(path: string): () => void {
  const fd = openSync(
    `${path}.lock`,
    constants.O_RDWR | constants.O_CREAT | constants.O_NOFOLLOW,
    0o600
  );
  try {
    lock_fd(fd, 'Checkpoint is in use by another process');
  } catch (error) {
    closeSync(fd);
    throw error;
  }
  // Keep the inode: unlinking it could let another writer acquire a different lock.
  return () => closeSync(fd);
}
export function write_all(fd: number, data: Buffer): void {
  let offset = 0;
  while (offset < data.length) {
    const count = writeSync(fd, data, offset, data.length - offset);
    if (!count) throw new Error('Output write made no progress');
    offset += count;
  }
}
export function atomic_json(path: string, data: unknown): void {
  const temporary = join(dirname(path), `${basename(path)}.${randomUUID()}.tmp`);
  const fd = openSync(temporary, 'wx', 0o600);
  try {
    try {
      write_all(fd, Buffer.from(`${json(data)}\n`));
      fsyncSync(fd);
    } finally {
      closeSync(fd);
    }
    renameSync(temporary, path);
    const directory = openSync(dirname(path), constants.O_RDONLY);
    try {
      fsyncSync(directory);
    } finally {
      closeSync(directory);
    }
  } finally {
    rmSync(temporary, { force: true });
  }
}
