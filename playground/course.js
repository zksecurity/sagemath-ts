import { createEditor } from './editor.js';
import { createRunner, resolveSage } from './runner.js';

const runtime = document.getElementById('runtime');
resolveSage(runtime);

for (const block of document.querySelectorAll('[data-example]')) {
  const textarea = block.querySelector('[data-editor]');
  const runButton = block.querySelector('[data-run]');
  const statusEl = block.querySelector('[data-status]');
  const outputEl = block.querySelector('[data-output]');

  if (!textarea || !runButton || !statusEl || !outputEl) continue;

  const editor = createEditor(textarea);
  const runner = createRunner({ outputEl, statusEl });

  // Get the cumulative prefix (all prior code blocks) from data attribute
  const prefix = block.dataset.prefix || '';

  runButton.addEventListener('click', () => {
    // Execute prefix + visible code for cumulative execution
    const fullCode = prefix + editor.getValue();
    runner.run(fullCode);
  });
}
