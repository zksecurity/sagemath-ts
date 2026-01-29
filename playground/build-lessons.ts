import { readdirSync, statSync, writeFileSync } from 'fs';
import { join } from 'path';
import { marked } from 'marked';
import { lessons, lessonNavMapping } from './lessons/lessons.generated.ts';

// Configure marked for safe rendering
marked.setOptions({
  breaks: true, // Convert \n to <br>
  gfm: true, // GitHub Flavored Markdown
});

const rootDir = new URL('..', import.meta.url).pathname;
const tutorialDir = join(rootDir, 'tutorial');
const outputPath = join(rootDir, 'playground', 'course.html');

// Use the nav mapping from generated lessons
const lessonAnchorMap = lessonNavMapping;

function escapeHtml(value: string): string {
  return value
    .replace(/&/g, '&amp;')
    .replace(/</g, '&lt;')
    .replace(/>/g, '&gt;')
    .replace(/"/g, '&quot;')
    .replace(/'/g, '&#39;');
}

function titleCaseFromSlug(slug: string): string {
  const small = new Set(['and', 'or', 'the', 'of', 'to', 'in', 'for', 'with']);
  return slug
    .split('-')
    .map((word, idx) => {
      if (idx > 0 && small.has(word)) return word;
      return word.charAt(0).toUpperCase() + word.slice(1);
    })
    .join(' ');
}

function parseChapterLabel(folder: string): string {
  const match = folder.match(/^(\d+)-(.*)$/);
  if (!match) return titleCaseFromSlug(folder);
  return `${match[1]}. ${titleCaseFromSlug(match[2])}`;
}

function parsePartTitle(folder: string): string {
  const match = folder.match(/^part(\d+)-(.*)$/);
  if (!match) return titleCaseFromSlug(folder);
  return `Part ${match[1]}: ${titleCaseFromSlug(match[2])}`;
}

function readTutorialNav() {
  const parts = readdirSync(tutorialDir)
    .filter((entry) => entry.startsWith('part') && statSync(join(tutorialDir, entry)).isDirectory())
    .sort();

  return parts.map((part) => {
    const partPath = join(tutorialDir, part);
    const chapters = readdirSync(partPath)
      .filter((entry) => statSync(join(partPath, entry)).isDirectory())
      .sort();

    return {
      title: parsePartTitle(part),
      items: chapters.map((chapter) => {
        const key = `${part}/${chapter}`;
        const anchor = lessonAnchorMap[key];
        return {
          label: parseChapterLabel(chapter),
          href: anchor ? `#lesson-${anchor}` : undefined,
          enabled: Boolean(anchor),
        };
      }),
    };
  });
}

const nav = readTutorialNav();

function renderNav() {
  return nav
    .map((group) => {
      const items = group.items
        .map((item) => {
          if (item.enabled && item.href) {
            return `<a href="${item.href}" class="nav-link">${escapeHtml(item.label)}</a>`;
          }
          return `<span class="nav-link muted">${escapeHtml(item.label)}</span>`;
        })
        .join('\n');
      return `\n<div class="nav-group">\n  <div class="nav-group-title">${escapeHtml(group.title)}</div>\n  ${items}\n</div>`;
    })
    .join('\n');
}

function renderMarkdown(text: string): string {
  // Render markdown to HTML
  const html = marked.parse(text) as string;
  // Wrap in a div with class for styling
  return `<div class="lesson-text">${html}</div>`;
}

function renderLessonBlocks(blocks: typeof lessons[number]['blocks']) {
  return blocks
    .map((block) => {
      if (block.type === 'paragraph') {
        return renderMarkdown(block.text);
      }
      if (block.type === 'code') {
        // Include prefix as a data attribute for cumulative execution
        const prefixAttr = block.prefix ? ` data-prefix="${escapeHtml(block.prefix)}"` : '';
        return `\n<div class="inline-example" data-example${prefixAttr}>\n  <div class="inline-header">\n    <span>Run it</span>\n    <button class="ghost" data-run>Run</button>\n    <span class="inline-status" data-status>Idle</span>\n  </div>\n  <textarea data-editor>${escapeHtml(block.code)}</textarea>\n  <pre data-output></pre>\n</div>`;
      }
      return '';
    })
    .join('\n');
}

function renderLessons() {
  return lessons
    .map((lesson) => {
      return `\n<section class="lesson" id="lesson-${lesson.id}">\n  <h2>${escapeHtml(lesson.title)}</h2>\n  ${renderLessonBlocks(lesson.blocks)}\n</section>`;
    })
    .join('\n');
}

const html = `<!doctype html>
<html lang="en">
  <head>
    <meta charset="UTF-8" />
    <meta name="viewport" content="width=device-width, initial-scale=1" />
    <title>sagemath-ts Lessons</title>
    <meta name="description" content="Textbook lessons with inline runnable examples." />
    <link rel="stylesheet" href="./styles.css" />
    <link rel="stylesheet" href="./course.css" />
  </head>
  <body>
    <div class="glow"></div>
    <div class="grain"></div>

    <nav class="top-nav">
      <div class="nav-links">
        <a href="./index.html">Home</a>
        <a href="./course.html" class="active">Lessons</a>
        <a href="./docs.html">Docs</a>
        <a href="./playground.html">Playground</a>
        <a href="https://github.com/ZkSecurity/sagemath-ts" target="_blank" rel="noreferrer">GitHub</a>
      </div>
    </nav>

    <header class="hero compact">
      <h1>Textbook lessons with runnable snippets.</h1>
      <p>
        Follow the narrative, tweak the code inline, and keep the playground
        nearby for experiments.
      </p>
      <div class="hero-actions">
        <div class="runtime" id="runtime">Runtime: loading...</div>
      </div>
    </header>

    <main class="course-layout">
      <aside class="course-sidebar">
        <div class="nav-title">Chapters</div>
        ${renderNav()}
      </aside>

      <article class="course">
        ${renderLessons()}
      </article>
    </main>

    <script type="module" src="./dist/sage-browser.js"></script>
    <script type="module" src="./dist/course.js"></script>
  </body>
</html>`;

writeFileSync(outputPath, html);
console.log(`Lessons written to ${outputPath}`);

// Validate code snippets run without throwing (warn only, don't fail build).
const Sage = await import('sagemath-ts');
const print = (..._args: unknown[]) => undefined;
let validationErrors = 0;
for (const lesson of lessons) {
  for (const block of lesson.blocks) {
    if (block.type !== 'code') continue;
    try {
      // Execute with prefix (cumulative code) + visible code
      const fullCode = (block.prefix || '') + block.code;
      const fn = new Function('Sage', 'print', 'console', fullCode);
      fn(Sage, print, console);
    } catch (err) {
      validationErrors++;
      console.warn(`⚠️  Lesson "${lesson.id}" has code that may not run: ${err instanceof Error ? err.message : err}`);
    }
  }
}
if (validationErrors > 0) {
  console.warn(`\n⚠️  ${validationErrors} code blocks have potential issues (missing functions, etc.)`);
}
