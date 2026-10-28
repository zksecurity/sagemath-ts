import { afterAll, beforeAll, expect, test } from 'bun:test';
import { spawnSync } from 'node:child_process';
import { generateKeyPairSync, randomInt } from 'node:crypto';
import { mkdtempSync, readFileSync, rmSync, writeFileSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import fc from 'fast-check';
import { candidate } from './core.js';
import { fingerprint, fingerprints, targets, valid_fingerprint } from './targets.js';

const dir = mkdtempSync(join(tmpdir(), 'timestamp-targets-'));
afterAll(() => rmSync(dir, { recursive: true, force: true }));
function openssl(args: string[]) {
  const result = spawnSync('openssl', args, { encoding: 'utf8' });
  if (result.status !== 0) throw new Error(result.stderr || String(result.error));
}
beforeAll(() => {
  for (const curve of ['prime256v1', 'secp256k1']) {
    openssl(['ecparam', '-name', curve, '-genkey', '-noout', '-out', join(dir, `${curve}.key`)]);
    openssl([
      'req',
      '-new',
      '-x509',
      '-key',
      join(dir, `${curve}.key`),
      '-out',
      join(dir, `${curve}.cert`),
      '-days',
      '1',
      '-subj',
      '/CN=synthetic test',
    ]);
    openssl([
      'x509',
      '-in',
      join(dir, `${curve}.cert`),
      '-outform',
      'DER',
      '-out',
      join(dir, `${curve}.der`),
    ]);
    openssl(['ec', '-in', join(dir, `${curve}.key`), '-pubout', '-out', join(dir, `${curve}.pub`)]);
    openssl([
      'ec',
      '-in',
      join(dir, `${curve}.key`),
      '-pubout',
      '-conv_form',
      'compressed',
      '-outform',
      'DER',
      '-out',
      join(dir, `${curve}.compressed`),
    ]);
    openssl([
      'pkcs8',
      '-topk8',
      '-nocrypt',
      '-in',
      join(dir, `${curve}.key`),
      '-outform',
      'DER',
      '-out',
      join(dir, `${curve}.private-der`),
    ]);
  }
});
function pem(der: Buffer): Buffer {
  return Buffer.from(
    `-----BEGIN PUBLIC KEY-----\n${der
      .toString('base64')
      .match(/.{1,64}/g)!
      .join('\n')}\n-----END PUBLIC KEY-----\n`
  );
}

test('live OpenSSL PEM/DER keys and certificates normalize across both supported curves', () => {
  for (const curve of ['prime256v1', 'secp256k1']) {
    const pub = readFileSync(join(dir, `${curve}.pub`));
    const expected = fingerprints(pub);
    for (const extension of ['cert', 'der', 'compressed'])
      expect(fingerprints(readFileSync(join(dir, `${curve}.${extension}`)))).toEqual(expected);
    expect(
      fingerprints(Buffer.concat([pub, readFileSync(join(dir, `${curve}.cert`)), pub]))
    ).toEqual(expected);
    expect(targets([join(dir, `${curve}.pub`), join(dir, `${curve}.der`)])).toEqual(expected);
  }
});

test('private, unsupported, malformed and trailing inputs are refused', () => {
  const { publicKey } = generateKeyPairSync('ed25519');
  const pub = readFileSync(join(dir, 'prime256v1.pub'));
  const der = readFileSync(join(dir, 'prime256v1.der'));
  for (const data of [
    readFileSync(join(dir, 'prime256v1.key')),
    readFileSync(join(dir, 'prime256v1.private-der')),
    publicKey.export({ format: 'pem', type: 'spki' }),
    Buffer.concat([der, Buffer.from('tail')]),
    Buffer.concat([pub, Buffer.from('garbage')]),
    Buffer.from('-----BEGIN PUBLIC KEY-----\nbroken'),
    Buffer.from('not a fingerprint'),
    Buffer.from([0xff]),
    Buffer.from([0x30, 0x80]),
    Buffer.from([0x30, 0x81, 0]),
  ])
    expect(() => fingerprints(Buffer.from(data))).toThrow();
  expect(() =>
    fingerprints(Buffer.concat([pub, readFileSync(join(dir, 'prime256v1.key'))]))
  ).toThrow('Only PUBLIC KEY');
  expect(() => fingerprint(generateKeyPairSync('ed25519').privateKey)).toThrow('public keys');
  const invalid = join(dir, 'invalid');
  writeFileSync(invalid, 'invalid');
  expect(() => targets(invalid)).toThrow('Invalid target file');
  expect(targets(null)).toBeNull();
  expect(valid_fingerprint('A'.repeat(64))).toBe(false);
});

test('fresh public points, compressed SPKI and text sets have identical fingerprints', () => {
  const replay = Number(process.env.WEAK_KEYS_SEED ?? randomInt(0, 0x7fffffff));
  console.info(`timestamp targets generator=v1 replay seed=${replay}`);
  fc.assert(
    fc.property(fc.bigInt({ min: 0n, max: (1n << 80n) - 1n }), (seed) => {
      for (const curve of ['p256', 'secp256k1'] as const) {
        const row = candidate(seed, 'seconds', 'cpp-mt19937', 'concat-be-mod', curve)!;
        const der = Buffer.from(row.spki_der_base64, 'base64');
        const compressed = Buffer.from(
          (curve === 'p256'
            ? '3039301306072a8648ce3d020106082a8648ce3d030107032200'
            : '3036301006072a8648ce3d020106052b8104000a032200') + row.sec1_compressed_hex,
          'hex'
        );
        const expected = new Set([row.spki_sha256]);
        expect(fingerprints(der)).toEqual(expected);
        expect(fingerprints(compressed)).toEqual(expected);
        expect(fingerprints(pem(compressed))).toEqual(expected);
        expect(
          fingerprints(
            Buffer.from(`# comment\n\n${row.spki_sha256.toUpperCase()}\n${row.spki_sha256}`)
          )
        ).toEqual(expected);
        expect(() => fingerprints(Buffer.concat([der, Buffer.from([0])]))).toThrow();
      }
    }),
    { seed: replay, numRuns: 20 }
  );
  expect(fingerprints(Buffer.from('# empty\n'))).toEqual(new Set());
});
