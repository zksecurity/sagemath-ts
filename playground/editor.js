import { EditorView, keymap } from '@codemirror/view';
import { EditorState } from '@codemirror/state';
import { javascript } from '@codemirror/lang-javascript';
import { defaultKeymap } from '@codemirror/commands';
import { lineNumbers } from '@codemirror/view';
import { oneDark } from '@codemirror/theme-one-dark';

export function createEditor(textarea) {
  let getValue = () => textarea.value;
  let setValue = (value) => {
    textarea.value = value;
  };
  let focus = () => textarea.focus();

  const api = {
    getValue: () => getValue(),
    setValue: (value) => setValue(value),
    focus: () => focus(),
  };

  const parent = textarea.parentElement;
  if (!parent) {
    return api;
  }

  const container = document.createElement('div');
  container.className = 'editor-cm';
  parent.insertBefore(container, textarea);

  try {
    const view = new EditorView({
      state: EditorState.create({
        doc: textarea.value,
        extensions: [
          lineNumbers(),
          keymap.of(defaultKeymap),
          javascript({ typescript: true }),
          oneDark, // includes theme + syntax highlighting
          EditorView.lineWrapping,
        ],
      }),
      parent: container,
    });

    textarea.style.display = 'none';
    container.style.display = 'block';

    getValue = () => view.state.doc.toString();
    setValue = (value) => {
      view.dispatch({
        changes: { from: 0, to: view.state.doc.length, insert: value },
      });
    };
    focus = () => view.focus();
  } catch (err) {
    console.error('CodeMirror init failed', err);
    container.remove();
  }

  return api;
}
