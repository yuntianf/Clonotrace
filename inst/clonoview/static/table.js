/* The patient and clone tables, in their own tab.
 *
 * WHY A SEPARATE PAGE. These are for putting beside the viewer and clicking through, so
 * a modal over the plots was the wrong shape: it hid the very thing the click changes.
 * As its own tab it can sit next to the viewer on a second half of the screen, and every
 * click updates the other tab in place.
 *
 * HOW THE TWO TABS TALK. A BroadcastChannel named "clonoview". This tab posts a
 * selection; the viewer applies it and posts its state back, which is what keeps the
 * highlighted row in step even when the selection was made in the viewer (by clicking a
 * point on the clone map, or stepping to the next clone). Same origin, so the access
 * cookie is already carried and no token handling is needed here.
 */

const $ = (id) => document.getElementById(id);
const fmt = (n) => n.toLocaleString('en-US');
const P = new URLSearchParams(location.search);

const T = {
  dataset: P.get('dataset') || null,
  view: P.get('view') === 'clones' ? 'clones' : 'patients',
  info: null, meta: null,
  rows: null,                 // clone rows, when loaded
  selLineage: null,           // what the viewer is currently showing
  selPatients: [],
  selLineages: [],            // several clones at once, unioned with selPatients
  /* Sort state per table. `key` names a column, `dir` is 1 or -1. Kept here rather than
   * in the DOM so it survives a re-render, which every selection click causes. */
  sort: { patients: { key: 'clones', dir: -1 }, clones: { key: 'cells', dir: -1 } },
};

/* ---- sorting -------------------------------------------------------------
 * Columns declare how to get their value; the comparator is shared. Numbers sort
 * numerically and text with localeCompare, decided from the value rather than declared,
 * because a column of patient ids like "02619-8" and "02619-40" sorts wrongly either as
 * pure text (8 after 40) or as a bare number (both NaN). `undefined` and null always sink
 * to the bottom whichever direction is chosen, so an empty cell never displaces real data
 * at the top of the table.
 */
function cmpVals(a, b) {
  const an = a === null || a === undefined || a === '', bn = b === null || b === undefined || b === '';
  if (an && bn) return 0;
  if (an) return 1;          // blanks sink, both directions
  if (bn) return -1;
  if (typeof a === 'number' && typeof b === 'number') return a - b;
  return String(a).localeCompare(String(b), 'en', { numeric: true, sensitivity: 'base' });
}

function sortRows(rows, cols, st) {
  const col = cols.find(c => c.key === st.key);
  if (!col) return rows;
  const out = [...rows];
  out.sort((x, y) => {
    const c = cmpVals(col.get(x), col.get(y));
    // Blanks are pinned to the bottom by cmpVals returning +/-1 regardless of order, so
    // the direction multiplier must not flip them back to the top.
    const xb = col.get(x) === null || col.get(x) === undefined || col.get(x) === '';
    const yb = col.get(y) === null || col.get(y) === undefined || col.get(y) === '';
    if (xb !== yb) return c;
    return c * st.dir;
  });
  return out;
}

/* A header cell that sorts. Clicking cycles descending -> ascending on the same column,
 * and starts descending on a new one (the useful default for counts). */
function sortableTh(label, key, st, rerender, numeric) {
  const th = el('th', numeric ? 'num' : null);
  const on = st.key === key;
  th.textContent = label + (on ? (st.dir < 0 ? ' \u25be' : ' \u25b4') : '');
  th.classList.add('sortable');
  if (on) th.setAttribute('aria-sort', st.dir < 0 ? 'descending' : 'ascending');
  th.title = `Sort by ${label}`;
  th.onclick = () => {
    if (st.key === key) st.dir = -st.dir;
    else { st.key = key; st.dir = -1; }
    rerender();
  };
  return th;
}



/* Same seam as app.js: offline, window.CV.route answers instead of the server. */
const LOCAL = () => (window.CV && window.CV.route) ? window.CV : null;

async function api(path) {
  const cv = LOCAL();
  if (cv) return cv.route(path);
  const r = await fetch(path);
  const j = await r.json().catch(() => ({ error: `HTTP ${r.status}` }));
  if (!r.ok || j.error) throw new Error(j.error || `HTTP ${r.status}`);
  return j;
}

let toastTimer = null;
function toast(msg, ms = 2600) {
  const t = $('toast');
  t.textContent = msg; t.hidden = false;
  clearTimeout(toastTimer);
  toastTimer = setTimeout(() => { t.hidden = true; }, ms);
}

function el(tag, cls, text) {
  const e = document.createElement(tag);
  if (cls) e.className = cls;
  if (text !== undefined) e.textContent = text;
  return e;
}

function send(msg) {
  if (!LINK.connected()) { toast('This browser cannot link the two tabs.'); return; }
  LINK.send(msg);
}

/* ---- patient table -------------------------------------------------------
 * Rows are patients, columns are timepoints. Each cell carries BOTH numbers -- clone
 * units and cells -- because neither implies the other and the run's cell threshold
 * sits between them: a patient can have thousands of cells at a visit and no clone that
 * reached ten of them.
 */
function renderPatients() {
  const body = $('t-body');
  body.innerHTML = '';
  $('t-hint').textContent =
    'clone units / cells at each timepoint. Click a row to show that patient\'s clones ' +
    'in the viewer — the left panel shows every clone unit with arrows between visits, ' +
    'the right shows their pooled cells and density. Shift-click to add patients to ' +
    'the selection; patients and individually chosen clones combine. Click any column ' +
    'heading to sort by it, again to reverse.';

  const visits = T.meta.visits;
  const fields = [...new Set(T.meta.patients.flatMap(p => Object.keys(p.fields)))].sort();
  const st = T.sort.patients;
  const cols = [
    { key: 'patient', label: 'Patient', get: p => p.id },
    ...fields.map(f => ({ key: 'f:' + f, label: f, get: p => p.fields[f] || null })),
    { key: 'clones', label: 'Clones', num: true, get: p => p.n_clone_units },
    { key: 'cells', label: 'Cells', num: true, get: p => p.n_cells },
    // Per-visit columns sort on clone units, the number the rest of the table leads
    // with; the cell count is shown beside it but is not what "most clones here" means.
    ...visits.map(v => ({ key: 'v:' + v.code, label: v.label, num: true,
                          get: p => p.per_visit[v.code].clone_units })),
  ];
  const t = el('table');
  const thead = el('thead'), hr = el('tr');
  for (const c of cols) hr.appendChild(sortableTh(c.label, c.key, st, renderPatients, c.num));
  thead.appendChild(hr); t.appendChild(thead);

  const tb = el('tbody');
  const rows = sortRows(T.meta.patients, cols, st);
  for (const p of rows) {
    const pi = T.info.patients.indexOf(p.id);
    const tr = el('tr');
    tr.setAttribute('aria-selected', T.selPatients.includes(pi) ? 'true' : 'false');
    tr.appendChild(el('td', null, p.id));
    for (const f of fields) tr.appendChild(el('td', null, p.fields[f] || '—'));
    tr.appendChild(el('td', 'num', fmt(p.n_clone_units)));
    tr.appendChild(el('td', 'num', fmt(p.n_cells)));
    for (const v of visits) {
      const pv = p.per_visit[v.code];
      const empty = pv.clone_units === 0 && pv.cells === 0;
      const td = el('td', 'num' + (empty ? ' zero' : ''));
      td.textContent = empty ? '—' : `${fmt(pv.clone_units)} / ${fmt(pv.cells)}`;
      // "Sampled, but nothing cleared the threshold" and "never sampled" are different
      // facts and the table must not let them look the same.
      if (pv.clone_units === 0 && pv.cells > 0) {
        td.classList.add('cell-count');
        td.title = `${fmt(pv.cells)} cells, but no clone reached the ` +
                   `${T.meta.expanded_min}-cell threshold at this visit`;
      }
      tr.appendChild(td);
    }
    tr.onclick = (ev) => {
      let sel;
      if (ev.shiftKey) {
        sel = T.selPatients.includes(pi)
          ? T.selPatients.filter(x => x !== pi) : [...T.selPatients, pi];
      } else {
        sel = T.selPatients.length === 1 && T.selPatients[0] === pi ? [] : [pi];
      }
      T.selPatients = sel;
      send((sel.length || T.selLineages.length)
        ? { type: 'select', patients: sel, lineages: T.selLineages }
        : { type: 'clear' });
      renderPatients();
    };
    tb.appendChild(tr);
  }
  t.appendChild(tb);
  body.appendChild(t);
}

/* ---- clone table ---------------------------------------------------------
 * One row per clone, stratified by patient. ONE click and the viewer changes.
 */
async function renderClones() {
  const body = $('t-body');
  body.innerHTML = 'loading…';
  const q = T.selPatients.length
    ? '?' + T.selPatients.map(p => `patient=${p}`).join('&') : '';
  let data;
  try { data = await api(`/api/lineages/${T.dataset}${q}`); }
  catch (e) { body.textContent = 'error: ' + e.message; return; }
  T.rows = data.rows;
  body.innerHTML = '';

  $('t-hint').innerHTML =
    `${fmt(data.n)} clones` +
    (T.selPatients.length
      ? ` for ${T.selPatients.map(p => T.info.patients[p]).join(', ')}` +
        ' — <b>use the Patients tab to change this</b>'
      : ' across all patients') +
    `. One click shows the clone in the viewer, shift-click adds it to the selection ` +
    `so several clones show at once (and they combine with any patients chosen on the ` +
    `Patients tab). Up to ${MARKERS.MAX} hand-picked clones each get their own marker ` +
    `shape on the clone map, shown here beside the row; colour still means profile. ` +
    `Click any column heading to sort, again to reverse. The squares are the sampling schedule ` +
    `(${T.meta.visits.map(v => v.code).join(' ')}): filled = a clone unit there, ` +
    `coloured by its profile; dashed = it had cells but fewer than ` +
    `${T.meta.expanded_min}, so Clonotrace gave it no position on the clone map; ` +
    `empty = no cells.`;

  const st = T.sort.clones;
  const cols = [
    { key: 'clone', label: 'Clone', get: r => r.label },
    { key: 'patient', label: 'Patient', get: r => T.info.patients[r.patient] },
    { key: 'span', label: 'Timepoints', num: true, get: r => r.visits.length },
    { key: 'units', label: 'Visits', num: true, get: r => r.n_units },
    { key: 'cells', label: 'Cells', num: true, get: r => r.n_cells },
    { key: 'profiles', label: 'Profiles', get: r => r.profiles[0] ?? null },
  ];
  const t = el('table');
  const thead = el('thead'), hr = el('tr');
  for (const c of cols)
    hr.appendChild(sortableTh(c.label, c.key, st, () => renderClones().catch(e => toast(e.message)), c.num));
  thead.appendChild(hr); t.appendChild(thead);

  const tb = el('tbody');
  for (const r of sortRows(data.rows, cols, st)) {
    const tr = el('tr');
    tr.dataset.lineage = r.lineage;
    tr.setAttribute('aria-selected',
      (T.selLineage === r.lineage || T.selLineages.includes(r.lineage)) ? 'true' : 'false');
    // The marker slot is filled by markSelected(), not here, so the shapes follow a
    // selection change without re-rendering and re-sorting the whole table.
    const idc = el('td');
    idc.appendChild(el('span', 'mk-slot'));
    idc.appendChild(el('span', null, r.label));
    tr.appendChild(idc);
    tr.appendChild(el('td', null, T.info.patients[r.patient]));

    const td = el('td'), mini = el('span', 'mini');
    for (let vi = 0; vi < T.meta.visits.length; vi++) {
      const i = el('i');
      const k = r.visits.indexOf(vi);
      if (k >= 0) {
        i.classList.add('on');
        i.style.background = T.info.profile_colours[r.profiles[k]];
        i.title = `${T.meta.visits[vi].label}: profile ${r.profiles[k]}`;
      } else if (r.visits_with_cells.includes(vi)) {
        i.classList.add('sub');
        i.title = `${T.meta.visits[vi].label}: cells present but below the ` +
                  `${T.meta.expanded_min}-cell threshold`;
      } else {
        i.title = `${T.meta.visits[vi].label}: no cells`;
      }
      mini.appendChild(i);
    }
    td.appendChild(mini); tr.appendChild(td);
    tr.appendChild(el('td', 'num', String(r.n_units)));
    tr.appendChild(el('td', 'num', fmt(r.n_cells)));
    tr.appendChild(el('td', null, r.profiles.join(' → ') || '—'));
    tr.onclick = (ev) => {
      if (ev.shiftKey) {
        // Accumulate. The viewer unions these with any selected patients, so a
        // "these three clones plus everything from patient 12" selection is reachable.
        T.selLineages = T.selLineages.includes(r.lineage)
          ? T.selLineages.filter(x => x !== r.lineage)
          : [...T.selLineages, r.lineage];
      } else {
        T.selLineages = [r.lineage];
        T.selPatients = [];
      }
      T.selLineage = T.selLineages.length === 1 ? T.selLineages[0] : null;
      send({ type: 'select', patients: T.selPatients, lineages: T.selLineages });
      markSelected();
    };
    tb.appendChild(tr);
  }
  t.appendChild(tb);
  body.appendChild(t);
  markSelected();
}

/* The shape the viewer is drawing this clone with, or null when it is not individually
 * marked. SORTED ASCENDING to match the server, which sorts `lineages_selected` before
 * the viewer indexes into it: our selection is in click order, and using that order
 * here would label the row with one shape while the map drew another. */
function markerForLineage(li) {
  const sel = [...T.selLineages].sort((a, b) => a - b);
  return MARKERS.shapeAt(sel.indexOf(li), sel.length);
}

function markSelected() {
  for (const tr of document.querySelectorAll('tbody tr[data-lineage]')) {
    const li = Number(tr.dataset.lineage);
    const on = li === T.selLineage || T.selLineages.includes(li);
    tr.setAttribute('aria-selected', on ? 'true' : 'false');
    const slot = tr.querySelector('.mk-slot');
    if (slot) {
      const shape = on ? markerForLineage(li) : null;
      slot.innerHTML = shape ? MARKERS.glyph(shape) : '';
      slot.title = shape ? 'This clone’s marker on the clone map' : '';
    }
    if (on && T.selLineages.length <= 1) tr.scrollIntoView({ block: 'nearest' });
  }
}

function setView(v) {
  T.view = v;
  $('tab-patients').setAttribute('aria-pressed', v === 'patients' ? 'true' : 'false');
  $('tab-clones').setAttribute('aria-pressed', v === 'clones' ? 'true' : 'false');
  if (v === 'patients') renderPatients();
  else renderClones().catch(e => toast(e.message));
}

async function boot() {
  if (!T.dataset) {
    const d = await api('/api/datasets');
    T.dataset = d.datasets.length ? d.datasets[0].name : null;
  }
  if (!T.dataset) { $('t-body').textContent = 'No prepared datasets.'; return; }
  T.info = await api(`/api/dataset/${T.dataset}`);
  T.meta = T.info.meta;
  $('t-dataset').textContent =
    `${T.meta.label} — ${fmt(T.meta.counts.patients)} patients, ` +
    `${fmt(T.meta.counts.lineages)} clones, ${fmt(T.meta.counts.clone_units)} clone units`;
  $('t-link').textContent = LINK.describe();

  $('tab-patients').onclick = () => setView('patients');
  $('tab-clones').onclick = () => setView('clones');
  $('t-clear').onclick = () => {
    T.selPatients = []; T.selLineage = null; T.selLineages = [];
    send({ type: 'clear' });
    setView(T.view);
  };

  LINK.onMessage((m) => {
    if (m.type !== 'state') return;
    // The viewer is the source of truth: it may have changed selection on its own,
    // by a click on the clone map or the next-clone button.
    const changedPatients =
      JSON.stringify(m.patients || []) !== JSON.stringify(T.selPatients);
    T.selLineage = m.lineage;
    T.selPatients = m.patients || [];
    T.selLineages = m.lineages || (m.lineage === null ? [] : [m.lineage]);
    if (T.view === 'clones') markSelected();
    if (changedPatients && T.view === 'patients') renderPatients();
  });
  LINK.send({ type: 'hello' });     // ask the viewer for its current state

  setView(T.view);
}

boot().catch(e => { $('t-body').textContent = 'error: ' + e.message; });
