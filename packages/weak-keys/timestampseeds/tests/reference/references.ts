/** Live native references. Missing compilers/runtimes fail; no saved-output fallback. */
import { execFileSync } from 'node:child_process';
import { mkdtempSync, rmSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { fileURLToPath } from 'node:url';
const ROOT = fileURLToPath(new URL('./', import.meta.url));
export class References {
  readonly directory = mkdtempSync(join(tmpdir(), 'timestamp-native-'));
  constructor(readonly windowsOnly = false) {
    try {
      for (const name of windowsOnly
        ? ['native_rng']
        : ['mt19937', 'lcg_examples', 'seeded_mt19937', 'native_rng']) {
        if (process.platform === 'win32')
          this.execute('cl', [
            '/nologo',
            '/EHsc',
            join(ROOT, `${name}.cpp`),
            `/Fe:${this.executable(name)}`,
            `/Fo:${join(this.directory, `${name}.obj`)}`,
          ]);
        else
          this.execute(process.env.CXX ?? 'c++', [
            '-std=c++11',
            join(ROOT, `${name}.cpp`),
            '-o',
            this.executable(name),
          ]);
      }
      if (!windowsOnly)
        this.execute(process.env.JAVAC ?? 'javac', [
          '-d',
          this.directory,
          join(ROOT, 'RandomReference.java'),
        ]);
    } catch (error) {
      this.close();
      throw error;
    }
  }
  executable(name: string): string {
    return join(this.directory, name + (process.platform === 'win32' ? '.exe' : ''));
  }
  execute(command: string, args: string[]): string {
    return execFileSync(command, args, {
      encoding: 'utf8',
      stdio: ['ignore', 'pipe', 'pipe'],
      maxBuffer: 8 * 1024 * 1024,
    });
  }
  java(seed: bigint, count: number, ...extra: string[]): string {
    return this.execute(process.env.JAVA ?? 'java', [
      '-cp',
      this.directory,
      'RandomReference',
      String(BigInt.asIntN(64, seed)),
      String(count),
      ...extra,
    ]);
  }
  words(name: string, seed: bigint, count: number): bigint[] {
    let text: string;
    if (name === 'java-random') text = this.java(seed, count);
    else if (name === 'python-mt19937')
      text = this.execute(this.executable('seeded_mt19937'), [String(seed), String(count)]);
    else if (name === 'cpp-mt19937')
      text = this.execute(this.executable('mt19937'), [
        String(BigInt.asUintN(32, seed)),
        String(count),
      ]);
    else {
      const executable =
        name === 'ansi-lcg31' || (name === 'msvc-rand' && !this.windowsOnly)
          ? 'lcg_examples'
          : 'native_rng';
      text = this.execute(this.executable(executable), [
        name,
        String(BigInt.asUintN(32, seed)),
        String(count),
      ]);
    }
    const values = text.trim().split(/\s+/).map(BigInt);
    if (values.length !== count) throw new Error('Native reference returned wrong word count');
    return values;
  }
  value(rng: string, seed: bigint, profile: string, skip: number, index: number): bigint {
    if (profile === 'java-nextbytes-reject') {
      const hex = this.java(seed, 32 * (index + 1), 'bytes', String(skip)).trim();
      return BigInt(`0x${hex.slice(-64)}`);
    }
    const width = rng === 'drand48-state' ? 6 : 4;
    const count = profile === 'low-byte-reject' ? 32 : Math.ceil(32 / width);
    const values = this.words(rng, seed, skip + (index + 1) * count).slice(skip + index * count);
    const little = [
      'concat-le-reject',
      'python-getrandbits-reject',
      'python-randbytes-reject',
    ].includes(profile);
    const material = Buffer.concat(
      values.map((value) => {
        if (profile === 'low-byte-reject') return Buffer.from([Number(value & 255n)]);
        const bytes = Buffer.from(value.toString(16).padStart(width * 2, '0'), 'hex');
        return little ? bytes.reverse() : bytes;
      })
    ).subarray(0, 32);
    if (profile === 'buffer-le-reject' || profile === 'python-getrandbits-reject')
      material.reverse();
    return BigInt(`0x${material.toString('hex')}`);
  }
  close(): void {
    rmSync(this.directory, { recursive: true, force: true });
  }
}
