function slugify(text) {
  return text
    .toLowerCase()
    .replace(/[^a-z0-9]+/g, '-')
    .replace(/(^-|-$)+/g, '');
}

function formatDate(iso) {
  const date = new Date(iso);
  if (Number.isNaN(date.getTime())) return iso;
  return date.toLocaleString();
}

function createFilterButtons(kinds, activeKind, onClick) {
  const wrapper = document.createElement('div');
  wrapper.className = 'docs-filter-buttons';
  for (const kind of kinds) {
    const button = document.createElement('button');
    button.type = 'button';
    button.dataset.kind = kind;
    button.textContent = kind === 'all' ? 'All' : kind;
    if (kind === activeKind) button.classList.add('active');
    button.addEventListener('click', () => onClick(kind));
    wrapper.appendChild(button);
  }
  return wrapper;
}

async function loadDocs() {
  const container = document.getElementById('docs-content');
  if (!container) return;

  try {
    const response = await fetch('./docs-data.json');
    if (!response.ok) throw new Error('docs-data.json not found');
    const data = await response.json();

    container.textContent = '';

    const layout = document.createElement('div');
    layout.className = 'docs-layout';

    const sidebar = document.createElement('aside');
    sidebar.className = 'docs-sidebar';

    const body = document.createElement('div');
    body.className = 'docs-body';

    const searchWrap = document.createElement('div');
    searchWrap.className = 'docs-search';
    const searchInput = document.createElement('input');
    searchInput.type = 'search';
    searchInput.placeholder = 'Search symbols, paths, or docs...';
    searchWrap.appendChild(searchInput);

    const docToggle = document.createElement('label');
    const docCheckbox = document.createElement('input');
    docCheckbox.type = 'checkbox';
    docCheckbox.checked = false;
    docToggle.appendChild(docCheckbox);
    docToggle.appendChild(document.createTextNode('Documented only'));

    const filterTitle = document.createElement('h3');
    filterTitle.textContent = 'Filters';

    const filters = document.createElement('div');
    filters.className = 'docs-filters';
    filters.appendChild(docToggle);

    const kinds = ['all'];
    for (const section of data.sections ?? []) {
      for (const item of section.items ?? []) {
        if (!kinds.includes(item.kind)) kinds.push(item.kind);
      }
    }

    let activeKind = 'all';
    const filterButtons = createFilterButtons(kinds, activeKind, (nextKind) => {
      activeKind = nextKind;
      for (const btn of filterButtons.querySelectorAll('button')) {
        btn.classList.toggle('active', btn.dataset.kind === activeKind);
      }
      applyFilters();
    });

    filters.appendChild(filterButtons);

    const meta = document.createElement('div');
    meta.className = 'docs-meta';
    meta.textContent = `Generated ${formatDate(data.generatedAt)} from ${data.source}.`;

    const navTitle = document.createElement('h3');
    navTitle.textContent = 'Modules';

    const nav = document.createElement('div');
    nav.className = 'docs-nav';

    const sectionElements = [];
    const navLinks = new Map();
    const itemElements = [];

    for (const section of data.sections ?? []) {
      const sectionId = `module-${slugify(section.title)}`;
      const link = document.createElement('a');
      link.href = `#${sectionId}`;
      link.textContent = section.title;
      link.dataset.section = sectionId;
      nav.appendChild(link);
      navLinks.set(sectionId, link);

      const sectionEl = document.createElement('section');
      sectionEl.className = 'docs-section';
      sectionEl.id = sectionId;

      const header = document.createElement('div');
      header.className = 'docs-section-header';
      const heading = document.createElement('h2');
      heading.textContent = section.title;
      header.appendChild(heading);
      const count = document.createElement('small');
      count.textContent = `${section.items.length} exports`;
      header.appendChild(count);

      sectionEl.appendChild(header);

      const itemsEl = document.createElement('div');
      itemsEl.className = 'docs-items';

      for (const item of section.items) {
        const card = document.createElement('article');
        card.className = 'docs-item';
        card.dataset.kind = item.kind;
        card.dataset.documented = item.doc ? 'true' : 'false';
        card.dataset.search = `${item.name} ${item.from} ${item.doc ?? ''}`.toLowerCase();

        const title = document.createElement('div');
        title.className = 'docs-item-title';
        const code = document.createElement('code');
        code.textContent = item.name;
        const kind = document.createElement('span');
        kind.className = 'docs-item-kind';
        kind.textContent = item.kind;
        title.appendChild(code);
        title.appendChild(kind);

        const metaLine = document.createElement('div');
        metaLine.className = 'docs-item-meta';
        metaLine.textContent = item.from;

        const doc = document.createElement('div');
        doc.className = 'docs-item-doc';
        if (item.doc) {
          doc.textContent = item.doc;
        } else {
          doc.classList.add('muted');
          doc.textContent = 'No docstring yet. Add JSDoc in the source to surface it here.';
        }

        card.appendChild(title);
        card.appendChild(metaLine);
        card.appendChild(doc);
        itemsEl.appendChild(card);
        itemElements.push(card);
      }

      sectionEl.appendChild(itemsEl);
      body.appendChild(sectionEl);
      sectionElements.push(sectionEl);
    }

    const emptyState = document.createElement('div');
    emptyState.className = 'docs-empty';
    emptyState.textContent = 'No matches. Try clearing filters.';
    emptyState.style.display = 'none';
    body.appendChild(emptyState);

    sidebar.appendChild(searchWrap);
    sidebar.appendChild(filterTitle);
    sidebar.appendChild(filters);
    sidebar.appendChild(meta);
    sidebar.appendChild(navTitle);
    sidebar.appendChild(nav);

    layout.appendChild(sidebar);
    layout.appendChild(body);
    container.appendChild(layout);

    function applyFilters() {
      const query = searchInput.value.trim().toLowerCase();
      const requireDocs = docCheckbox.checked;

      let visibleItems = 0;
      for (const item of itemElements) {
        const matchesQuery = !query || item.dataset.search?.includes(query);
        const matchesDocs = !requireDocs || item.dataset.documented === 'true';
        const matchesKind = activeKind === 'all' || item.dataset.kind === activeKind;
        const visible = matchesQuery && matchesDocs && matchesKind;
        item.style.display = visible ? '' : 'none';
        if (visible) visibleItems += 1;
      }

      let visibleSections = 0;
      for (const section of sectionElements) {
        const hasVisible = Array.from(section.querySelectorAll('.docs-item')).some(
          (item) => item.style.display !== 'none',
        );
        section.style.display = hasVisible ? '' : 'none';
        if (hasVisible) {
          visibleSections += 1;
          navLinks.get(section.id)?.removeAttribute('aria-hidden');
          navLinks.get(section.id)?.style.removeProperty('display');
        } else {
          navLinks.get(section.id)?.setAttribute('aria-hidden', 'true');
          navLinks.get(section.id)?.style.setProperty('display', 'none');
        }
      }

      emptyState.style.display = visibleItems === 0 ? 'block' : 'none';
      meta.textContent = `Showing ${visibleItems} exports across ${visibleSections} modules. Generated ${formatDate(
        data.generatedAt,
      )} from ${data.source}.`;
    }

    searchInput.addEventListener('input', applyFilters);
    docCheckbox.addEventListener('change', applyFilters);

    applyFilters();

    if ('IntersectionObserver' in window) {
      const observer = new IntersectionObserver(
        (entries) => {
          for (const entry of entries) {
            if (!entry.isIntersecting) continue;
            const id = entry.target.id;
            for (const link of navLinks.values()) {
              link.classList.toggle('active', link.dataset.section === id);
            }
          }
        },
        { rootMargin: '-40% 0px -50% 0px', threshold: 0.1 },
      );

      for (const section of sectionElements) {
        observer.observe(section);
      }
    }
  } catch (err) {
    container.textContent = 'Unable to load docs.';
  }
}

loadDocs();
