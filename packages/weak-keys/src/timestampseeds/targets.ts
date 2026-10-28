/** Normalize local public inputs to SHA256 of canonical, uncompressed SPKI DER. */
import { type KeyObject, X509Certificate, createHash } from 'node:crypto';
import { readFileSync } from 'node:fs';
import { validateEcKeys } from '../../../zksecurity-cheatsheets/src/attacks/ec-keys.js';
import { p256, secp256k1 } from '../../../zksecurity-cheatsheets/src/curves/params.js';
import { CURVES } from './core.js';

export function valid_fingerprint(value: unknown): value is string {
  return typeof value === 'string' && /^[0-9a-f]{64}$/.test(value);
}
interface Der {
  tag: number;
  body: Buffer;
  raw: Buffer;
}
/** Definite, minimal DER lengths with exact boundaries, including nested containers. */
function elements(data: Buffer): Der[] {
  const result: Der[] = [];
  let offset = 0;
  while (offset < data.length) {
    const start = offset;
    if (data.length - offset < 2) throw new Error('Truncated DER input');
    const tag = data[offset++]!;
    let length = data[offset++]!;
    if (length & 128) {
      const count = length & 127;
      if (!count || count > 4 || offset + count > data.length || data[offset] === 0)
        throw new Error('Invalid DER length');
      length = data.readUIntBE(offset, count);
      offset += count;
      if (length < 128) throw new Error('Invalid DER length');
    }
    if (offset + length > data.length) throw new Error('Truncated DER input');
    result.push({
      tag,
      body: data.subarray(offset, offset + length),
      raw: data.subarray(start, offset + length),
    });
    offset += length;
  }
  return result;
}
function sequence(data: Buffer): Der[] {
  const outer = elements(data);
  if (outer.length !== 1 || outer[0]!.tag !== 0x30)
    throw new Error('DER input must contain exactly one object');
  return elements(outer[0]!.body);
}
export function fingerprint(publicKey: KeyObject | Buffer): string {
  if (!Buffer.isBuffer(publicKey)) {
    if (publicKey.type !== 'public')
      throw new Error('Targets must use P-256 or secp256k1 public keys');
    return fingerprint(publicKey.export({ format: 'der', type: 'spki' }));
  }
  const fields = sequence(publicKey);
  if (
    fields.length !== 2 ||
    fields[0]!.tag !== 0x30 ||
    fields[1]!.tag !== 3 ||
    fields[1]!.body[0] !== 0
  )
    throw new Error('Invalid public key SPKI');
  const algorithm = elements(fields[0]!.body);
  if (
    algorithm.length !== 2 ||
    algorithm.some((item) => item.tag !== 6) ||
    algorithm[0]!.body.toString('hex') !== '2a8648ce3d0201'
  )
    throw new Error('Targets must use P-256 or secp256k1 public keys');
  const oid = algorithm[1]!.body.toString('hex');
  const curve = oid === '2a8648ce3d030107' ? 'p256' : oid === '2b8104000a' ? 'secp256k1' : null;
  if (!curve) throw new Error('Targets must use P-256 or secp256k1 public keys');
  const bytes = fields[1]!.body.subarray(1);
  if (
    !(
      (bytes.length === 33 && (bytes[0] === 2 || bytes[0] === 3)) ||
      (bytes.length === 65 && bytes[0] === 4)
    )
  )
    throw new Error('Invalid SEC1 public point');
  const checked = validateEcKeys(curve === 'p256' ? p256 : secp256k1, [bytes])[0]!;
  if (!checked.valid || !checked.point) throw new Error('Invalid public point');
  const point = checked.point;
  const der = Buffer.from(
    CURVES[curve].prefix +
      '04' +
      point.x.toString(16).padStart(64, '0') +
      point.y.toString(16).padStart(64, '0'),
    'hex'
  );
  return createHash('sha256').update(der).digest('hex');
}
function certificate(data: Buffer): string {
  // Validate the certificate container without asking the host crypto library to
  // decode its EC key (Bun's BoringSSL does not support secp256k1).
  new X509Certificate(data);
  const cert = sequence(data);
  if (cert.length !== 3 || cert[0]!.tag !== 0x30 || cert[1]!.tag !== 0x30 || cert[2]!.tag !== 3)
    throw new Error('Invalid certificate');
  const tbs = elements(cert[0]!.body);
  const spki = tbs[tbs[0]?.tag === 0xa0 ? 6 : 5];
  if (!spki) throw new Error('Invalid certificate');
  return fingerprint(spki.raw);
}
function pem_der(text: string): Buffer {
  const encoded = text.replace(/-----[^\n]+-----/g, '').replace(/\s/g, '');
  if (!/^(?:[A-Za-z0-9+/]{4})*(?:[A-Za-z0-9+/]{2}==|[A-Za-z0-9+/]{3}=)?$/.test(encoded))
    throw new Error('Malformed PEM input');
  return Buffer.from(encoded, 'base64');
}
/** @see Deviation: Weak-key application package (no SageMath counterpart) */
export function fingerprints(data: Buffer): Set<string> {
  if (data.includes('-----BEGIN ')) {
    const text = data.toString('utf8');
    const result = new Set<string>();
    let end = 0;
    for (const match of text.matchAll(/-----BEGIN ([A-Z0-9 ]+)-----\s+[\s\S]*?-----END \1-----/g)) {
      if (text.slice(end, match.index).trim())
        throw new Error('Unexpected content outside PEM blocks');
      const label = match[1];
      if (label === 'PUBLIC KEY') result.add(fingerprint(pem_der(match[0])));
      else if (label === 'CERTIFICATE') result.add(certificate(pem_der(match[0])));
      else throw new Error('Only PUBLIC KEY and CERTIFICATE PEM blocks are accepted');
      end = match.index! + match[0].length;
    }
    if (!result.size || text.slice(end).trim()) throw new Error('Malformed PEM input');
    return result;
  }
  const textLike = data.every((byte) => [9, 10, 13].includes(byte) || (byte >= 32 && byte < 127));
  if (data[0] === 0x30 && !textLike) {
    sequence(data);
    try {
      return new Set([fingerprint(data)]);
    } catch {
      return new Set([certificate(data)]);
    }
  }
  if (data.some((byte) => byte > 127)) throw new Error('Expected ASCII fingerprints');
  const result = new Set<string>();
  let text = data.toString('ascii');
  for (const code of [0x1c, 0x1d, 0x1e]) text = text.replaceAll(String.fromCharCode(code), '\n');
  for (let line of text.split(/\r\n|[\n\r\v\f]/)) {
    line = line.trim().toLowerCase();
    if (!line || line.startsWith('#')) continue;
    if (!valid_fingerprint(line))
      throw new Error(
        'Expected SHA256 SPKI fingerprints, PEM public keys/certificates, or a DER public key/certificate'
      );
    result.add(line);
  }
  return result;
}
export function targets(paths?: string | string[] | null): Set<string> | null {
  if (paths == null) return null;
  const result = new Set<string>();
  for (const path of typeof paths === 'string' ? [paths] : paths) {
    try {
      for (const hash of fingerprints(readFileSync(path))) result.add(hash);
    } catch (error) {
      throw new Error(`Invalid target file ${path}: ${(error as Error).message}`);
    }
  }
  return result;
}
