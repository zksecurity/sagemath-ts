// Synthetic test-only reference. Requires Node/OpenSSL, not Bun's crypto compatibility layer.
const { createECDH } = require('node:crypto');
const curve = { p256: 'prime256v1', p384: 'secp384r1', secp256k1: 'secp256k1' }[process.argv[2]];
const scalar = Number(process.argv[3]);
if (!curve || !Number.isSafeInteger(scalar) || scalar <= 0 || scalar > 65536) process.exit(2);
const key = createECDH(curve);
const secret = Buffer.alloc(32);
secret.writeUInt32BE(scalar, 28);
key.setPrivateKey(secret);
process.stdout.write(key.getPublicKey('hex', 'uncompressed'));
