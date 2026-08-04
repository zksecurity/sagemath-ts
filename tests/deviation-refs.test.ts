/**
 * Guards `@see Deviation:` tags against silently drifting away from `DEVIATIONS.md`.
 *
 * The tag is used two ways in this codebase: as a *reference* to a `DEVIATIONS.md`
 * section, and as inline prose describing the divergence at the call site. Both are
 * fine. What is not fine is a reference to a section that was later renamed — the tag
 * still reads like a working pointer, but leads nowhere.
 *
 * So this test does not demand that every tag name a section. It flags two shapes that
 * can only be broken references:
 *
 *   1. A tag that is a *near* match for a heading without being an exact one — what a
 *      section rename leaves behind. Measured on this tree, real renames score ~0.71 and
 *      the closest prose tag scores 0.53, so the threshold sits at 0.65.
 *   2. A kebab-case slug (`no-number-coercion`), which is only ever written as a pointer
 *      and scores too low to be caught by (1).
 */

import { describe, expect, test } from 'bun:test';
import { readFileSync, readdirSync, statSync } from 'node:fs';
import { join } from 'node:path';

const projectRoot = new URL('..', import.meta.url).pathname;

/** Similarity in [0, 1]: longest-common-subsequence length over the longer string. */
function similarity(a: string, b: string): number {
  const s = a.toLowerCase();
  const t = b.toLowerCase();
  if (s === t) return 1;
  if (s.length === 0 || t.length === 0) return 0;
  let prev = new Array<number>(t.length + 1).fill(0);
  for (let i = 1; i <= s.length; i++) {
    const cur = new Array<number>(t.length + 1).fill(0);
    for (let j = 1; j <= t.length; j++) {
      cur[j] = s[i - 1] === t[j - 1] ? prev[j - 1]! + 1 : Math.max(prev[j]!, cur[j - 1]!);
    }
    prev = cur;
  }
  return prev[t.length]! / Math.max(s.length, t.length);
}

function collectTsFiles(dir: string, out: string[] = []): string[] {
  for (const entry of readdirSync(dir)) {
    if (entry === 'node_modules' || entry === 'dist' || entry.startsWith('.')) continue;
    const full = join(dir, entry);
    if (statSync(full).isDirectory()) collectTsFiles(full, out);
    else if (entry.endsWith('.ts')) out.push(full);
  }
  return out;
}

describe('DEVIATIONS.md cross-references', () => {
  const deviations = readFileSync(join(projectRoot, 'DEVIATIONS.md'), 'utf8');
  const headings = deviations
    .split('\n')
    .map((line) => /^#{2,4}\s+(.*?)\s*$/.exec(line)?.[1])
    .filter((h): h is string => Boolean(h));

  test('DEVIATIONS.md has headings to reference', () => {
    expect(headings.length).toBeGreaterThan(50);
  });

  test('no @see Deviation tag is a dangling reference', () => {
    const files = collectTsFiles(join(projectRoot, 'packages'));
    const broken: string[] = [];

    for (const file of files) {
      const lines = readFileSync(file, 'utf8').split('\n');
      lines.forEach((line, i) => {
        const tag = /@see Deviation:\s*(.+?)\s*$/
          .exec(line)?.[1]
          ?.replace(/[*/\s]+$/, '')
          .trim();
        if (!tag) return;
        // An exact hit, or text lifted verbatim from the document, is fine.
        if (headings.some((h) => h.toLowerCase() === tag.toLowerCase())) return;
        if (deviations.includes(tag)) return;

        const rel = file.slice(projectRoot.length);

        // A kebab-case slug is always a pointer, never prose.
        if (/^[a-z0-9]+(-[a-z0-9]+)+$/.test(tag)) {
          broken.push(`${rel}:${i + 1}\n    slug "${tag}" matches no DEVIATIONS.md heading`);
          return;
        }

        const best = headings.reduce(
          (acc, h) => {
            const score = similarity(tag, h);
            return score > acc.score ? { score, heading: h } : acc;
          },
          { score: 0, heading: '' }
        );
        if (best.score >= 0.65) {
          broken.push(`${rel}:${i + 1}\n    tag:     "${tag}"\n    heading: "${best.heading}"`);
        }
      });
    }

    expect(broken).toEqual([]);
  });
});
