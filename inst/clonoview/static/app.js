/* ClonoView front end.
 *
 * THE ONE GEOMETRIC IDEA. Both panels draw into the same fixed data extent that the
 * bundle stores, and both fit it the same way: preserve aspect, centre, letterbox. The
 * cell panel's backdrop is a server-rendered PNG whose figure was sized to exactly that
 * extent's aspect, so an <svg> with viewBox = extent and the default
 * preserveAspectRatio="xMidYMid meet" lands pixel-for-pixel on top of an <img> with
 * object-fit: contain. The clone panel is a canvas, so it applies the same fit by hand
 * in fitBox(). Nothing has to be negotiated between the two sides -- if a future change
 * stops using meta.extent, every overlay silently shifts, which is why both panels take
 * it from the same place.
 *
 * SVG y grows downward and the embedding's y grows upward, so every drawn point is
 * (x, -y) and the viewBox is anchored at -y1. That single sign flip is the whole of it.
 */

const $ = (id) => document.getElementById(id);

const S = {
  dataset: null, meta: null, info: null,
  clones: null,            // clone-unit arrays
  cloneColour: true, cloneLabels: true,
  layout: 'default',
  axis: 'time', fixed: -1,
  /* One selection, two kinds of member. Patients and individual clones are
   * unioned by /api/set, so they can be held at the same time; a lone clone
   * keeps the richer single-clone view. */
  patientsSel: [], lineagesSel: [],
  cellMode: 'celltype', cellLabels: true,
  // Two independent selections, on purpose:
  //   profileSel  filters WHICH cells and clones are shown, in both panels
  //   isolate     highlights one class of whatever the cell panel is coloured by
  // Collapsing them into one control would make "profile 7 coloured by cell type"
  // impossible to ask for, and that is the view that says what a profile consists of.
  profileSel: null,
  isolate: null,           // {kind:'profile'|'celltype', value:int} or null
  patientSel: null,        // index into info.patients, or null for all
  clone: null,             // the selected lineage, with all its frames
  set: null,               // a selected SET of clones (one or more patients)
  frame: 0,                // which frame of S.clone is showing
  lineageRows: null,       // the clone table's current rows, for "next clone"
  mana: false,             // MANAscore overlay, on both panels at once
};

/* ---- tiny helpers -------------------------------------------------------- */

/* THE ONE SEAM BETWEEN SERVED AND OFFLINE. When this page is opened from file:// there is
 * no server to answer, and fetch is refused across the opaque origin a file:// page gets;
 * app/api.js then defines window.CV.route, which computes the same objects from the
 * packed arrays in the browser. Everything below this function is identical in both
 * modes, deliberately: the two versions must not become two front ends, because the
 * served one is what the offline one is checked against. */
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
function toast(msg, ms = 3200) {
  const t = $('toast');
  t.textContent = msg; t.hidden = false;
  clearTimeout(toastTimer);
  toastTimer = setTimeout(() => { t.hidden = true; }, ms);
}

const fmt = (n) => n.toLocaleString('en-US');

/* Fit a data extent into a pixel box, preserving aspect and centring -- the manual
 * equivalent of object-fit: contain / preserveAspectRatio="xMidYMid meet". */
function fitBox(extent, w, h) {
  const [x0, x1, y0, y1] = extent;
  const s = Math.min(w / (x1 - x0), h / (y1 - y0));
  const ox = (w - (x1 - x0) * s) / 2, oy = (h - (y1 - y0) * s) / 2;
  return {
    s, ox, oy,
    X: (x) => ox + (x - x0) * s,
    Y: (y) => h - oy - (y - y0) * s,       // flip: data y up, canvas y down
    invX: (px) => x0 + (px - ox) / s,
    invY: (py) => y0 + (h - oy - py) / s,
  };
}

/* ---- clone panel (canvas, ~1.5k points, fully interactive) --------------- */
function drawClones() {
  const cv = $('clone-canvas'), plot = $('clone-plot');
  if (!S.clones) return;
  const w = plot.clientWidth, h = plot.clientHeight;
  const dpr = window.devicePixelRatio || 1;
  cv.width = Math.round(w * dpr); cv.height = Math.round(h * dpr);
  const g = cv.getContext('2d');
  g.setTransform(dpr, 0, 0, dpr, 0, 0);
  g.clearRect(0, 0, w, h);
  g.fillStyle = S.info.colours.surface; g.fillRect(0, 0, w, h);

  const C = S.clones, f = fitBox(C.extent, w, h);
  S.cloneFit = f;
  const cols = S.info.profile_colours;
  const ctx = S.info.colours.context;
  const iso = S.profileSel;

  // Context first, highlighted on top, so a selected clone is never buried. A point is
  // "on" only if it survives every active filter -- profile, patient, and (when a clone
  // is open) not being that clone, whose own marks are drawn separately and larger.
  const r = 2.6;
  const pass = [[], []];
  for (let i = 0; i < C.x.length; i++) {
    let on = iso === null || C.profile[i] === iso;
    if (on && S.patientSel !== null) on = C.patient[i] === S.patientSel;
    if (on && active()) on = false;     // the whole map is context behind a selection
    pass[on ? 1 : 0].push(i);
  }
  for (const i of pass[0]) {
    g.fillStyle = ctx;
    g.beginPath(); g.arc(f.X(C.x[i]), f.Y(C.y[i]), r, 0, 7); g.fill();
  }
  for (const i of pass[1]) {
    g.fillStyle = S.cloneColour ? (cols[C.profile[i]] || ctx) : '#6d6b64';
    g.beginPath(); g.arc(f.X(C.x[i]), f.Y(C.y[i]), r, 0, 7); g.fill();
  }

  drawPtrcRings(g, f);
  if (active()) drawTrack(g, f);

  drawCloneLabels(f);
  if (!active()) {
    const bits = [`${fmt(C.x.length)} clone units`];
    if (S.patientSel !== null) {
      bits.push(`${S.info.patients[S.patientSel]}: ` +
                `${fmt(pass[1].length)} shown`);
    }
    if (iso !== null) bits.push(`profile ${iso} isolated`);
    $('clone-sub').textContent = bits.join(' — ');
  }
}

/* MANAscore-positive clones, ringed in red.
 *
 * Drawn as a THICK RING rather than a fill, for two reasons. It has to survive on top of
 * whatever the point is already coloured by -- profile colour carries the Clonotrace
 * answer and must not be overwritten by the MANAscore one, since the whole point of
 * showing them together is to see whether two independent methods agree. And a ring is a
 * different visual channel from hue, so it still reads where the categorical palette
 * already contains a red.
 *
 * A ring means pTRC = 1. No ring means either called-and-negative or, for the clones
 * whose patient had no usable cutoff, not callable at all -- those are drawn with a
 * dashed grey ring instead, because "not called" is not "called negative".
 */
function drawPtrcRings(g, f) {
  if (!S.mana || !S.clones || !S.clones.ptrc) return;
  const C = S.clones;
  const red = (S.info.colours && S.info.colours.mana) || '#e00025';
  for (let i = 0; i < C.x.length; i++) {
    if (S.patientSel !== null && C.patient[i] !== S.patientSel) continue;
    if (S.profileSel !== null && C.profile[i] !== S.profileSel) continue;
    const px = f.X(C.x[i]), py = f.Y(C.y[i]);
    if (C.ptrc[i] === 1) {
      g.save();
      g.beginPath(); g.arc(px, py, 6.0, 0, 7);
      g.strokeStyle = red; g.lineWidth = 2.6; g.stroke();
      g.restore();
    } else if (C.mana_usable && C.mana_usable[i] === 0) {
      g.save();
      g.beginPath(); g.arc(px, py, 5.2, 0, 7);
      g.strokeStyle = '#898781'; g.lineWidth = 1.0; g.setLineDash([2, 2]); g.stroke();
      g.restore();
    }
  }
}

/* Direct labels at each profile's centroid. Mandatory relief, not decoration: four of
 * the eleven profile colours sit below 3:1 contrast on this surface, so identity must
 * not rest on colour alone. Median, not mean -- a diffuse profile's mean lands in a
 * hole where none of its clones are. */
function drawCloneLabels(f) {
  const svg = $('clone-overlay'), plot = $('clone-plot');
  const w = plot.clientWidth, h = plot.clientHeight;
  svg.setAttribute('viewBox', `0 0 ${w} ${h}`);
  svg.innerHTML = '';
  if (!S.cloneLabels || !S.clones) return;
  const C = S.clones, cols = S.info.profile_colours;
  const by = new Map();
  for (let i = 0; i < C.x.length; i++) {
    if (!by.has(C.profile[i])) by.set(C.profile[i], []);
    by.get(C.profile[i]).push(i);
  }
  const iso = S.profileSel;
  for (const [p, idx] of [...by.entries()].sort((a, b) => a[0] - b[0])) {
    if (iso !== null && p !== iso) continue;
    const xs = idx.map(i => C.x[i]).sort((a, b) => a - b);
    const ys = idx.map(i => C.y[i]).sort((a, b) => a - b);
    const mx = xs[xs.length >> 1], my = ys[ys.length >> 1];
    svg.appendChild(chip(f.X(mx), f.Y(my), String(p), cols[p]));
  }
}

function chip(px, py, text, stroke) {
  const g = document.createElementNS('http://www.w3.org/2000/svg', 'g');
  const wid = 9 + text.length * 6.2;
  const rect = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
  rect.setAttribute('x', px - wid / 2); rect.setAttribute('y', py - 8);
  rect.setAttribute('width', wid); rect.setAttribute('height', 16);
  rect.setAttribute('rx', 5);
  rect.setAttribute('fill', 'rgba(252,252,251,.93)');
  rect.setAttribute('stroke', stroke || '#0b0b0b');
  rect.setAttribute('stroke-width', '1');
  const t = document.createElementNS('http://www.w3.org/2000/svg', 'text');
  t.setAttribute('x', px); t.setAttribute('y', py + 4);
  t.setAttribute('text-anchor', 'middle');
  t.setAttribute('font-size', '10.5'); t.setAttribute('font-weight', '650');
  t.setAttribute('fill', '#0b0b0b');
  t.textContent = text;
  g.appendChild(rect); g.appendChild(t);
  return g;
}

/* ---- cell panel (server raster + client overlay) ------------------------- */
let bgToken = 0;
function updateCellBg() {
  const img = $('cell-bg'), plot = $('cell-plot');
  const width = Math.min(2000, Math.max(700,
    Math.round(plot.clientWidth * (window.devicePixelRatio || 1))));
  const iso = S.isolate;
  const sel = (iso && iso.kind === (S.cellMode === 'profile' ? 'profile' : 'celltype'))
    ? iso.value : null;
  let q = `mode=${S.cellMode}&width=${width}`;
  if (S.profileSel !== null) q += `&profile=${S.profileSel}`;
  if (S.mana) q += '&mana=1';
  if (sel !== null) q += `&select=${sel}`;
  const url = `/api/bg/${S.dataset}?${q}`;
  if (img.dataset.url === url) { drawCellLabels(); return; }
  const my = ++bgToken;
  $('cell-spin').hidden = false;

  /* Offline, the same picture is drawn here instead of fetched. The URL is still built
   * above and still used as the cache key, so the "did anything actually change" test is
   * one piece of code rather than two -- and the <img> stays an <img>, which is what
   * keeps the SVG overlay's pixel-for-pixel alignment (object-fit: contain against
   * preserveAspectRatio) working without a second geometry path to get wrong. */
  const cv = LOCAL();
  if (cv) {
    const b = cv.loadedDataset(S.dataset);
    if (!b) { $('cell-spin').hidden = true; return; }
    // Yielded to the event loop first: rasterising 273k cells is tens of milliseconds of
    // straight-line work, and running it inside the click handler freezes the control
    // that started it, which reads as the app having hung rather than being busy.
    setTimeout(() => {
      if (my !== bgToken) return;
      let canvas;
      try {
        canvas = cv.background(b, {
          mode: S.cellMode, width, select: sel,
          profile: S.profileSel, mana: S.mana,
        });
      } catch (e) {
        $('cell-spin').hidden = true;
        toast(`could not draw the cell embedding: ${e.message}`);
        return;
      }
      canvas.toBlob((blob) => {
        if (my !== bgToken || !blob) return;   // a newer redraw won; nothing allocated yet
        const old = img.dataset.blob;
        const next = URL.createObjectURL(blob);
        img.src = next; img.dataset.url = url; img.dataset.blob = next;
        // Object URLs are held until revoked; one per redraw would leak a megabyte a
        // click over a session of stepping through filters.
        if (old) URL.revokeObjectURL(old);
        $('cell-spin').hidden = true;
        drawCellLabels();
      }, 'image/png');
    }, 0);
    return;
  }

  const pre = new Image();
  pre.onload = () => {
    if (my !== bgToken) return;              // a newer request won
    img.src = url; img.dataset.url = url;
    $('cell-spin').hidden = true;
    drawCellLabels();
  };
  pre.onerror = () => {
    if (my !== bgToken) return;
    $('cell-spin').hidden = true;
    toast('could not render the cell embedding');
  };
  pre.src = url;
}

function drawCellLabels() {
  const svg = $('cell-overlay');
  const ext = S.meta.extent.cell;
  svg.setAttribute('viewBox', `${ext[0]} ${-ext[3]} ${ext[1] - ext[0]} ${ext[3] - ext[2]}`);
  svg.innerHTML = '';
  if (!active()) $('cell-sub').textContent = cellSubtitle();
  // The clone overlay shares this <svg>, which is cleared above -- so it has to be
  // redrawn here rather than called separately, or an async background load would wipe
  // the clone off the map while leaving the timebar saying it is shown.
  const drawCloneAfter = () => { drawCellSelection(); };
  if (!S.cellLabels || S.cellMode === 'plain') { drawCloneAfter(); return; }

  // Centroids come from S.centroids, which is refetched whenever the profile filter
  // changes -- a label placed at the unfiltered centroid can land where the selected
  // profile has no cells at all.
  const src = S.centroids || S.info.centroids;
  const cents = S.cellMode === 'profile' ? src.profile : src.cell_type;
  const cols = S.cellMode === 'profile'
    ? S.info.profile_colours : S.info.cell_type_colours;
  const names = S.cellMode === 'profile' ? null : S.info.cell_types;
  const iso = S.isolate && S.isolate.kind ===
    (S.cellMode === 'profile' ? 'profile' : 'celltype') ? S.isolate.value : null;

  // Labels are drawn in data units here, so they scale with the viewBox. Sized from the
  // extent rather than in pixels for that reason.
  const span = ext[1] - ext[0];
  const fs = span / 46;
  for (const [k, xy] of Object.entries(cents)) {
    const key = Number(k);
    if (iso !== null && key !== iso) continue;
    const label = names ? names[key] : String(key);
    const t = document.createElementNS('http://www.w3.org/2000/svg', 'text');
    t.setAttribute('x', xy[0]); t.setAttribute('y', -xy[1]);
    t.setAttribute('text-anchor', 'middle');
    t.setAttribute('font-size', fs);
    t.setAttribute('font-weight', '650');
    t.setAttribute('fill', '#0b0b0b');
    t.setAttribute('stroke', 'rgba(252,252,251,.92)');
    t.setAttribute('stroke-width', fs * 0.28);
    t.setAttribute('paint-order', 'stroke');
    t.textContent = label;
    svg.appendChild(t);
  }
  drawCloneAfter();
}

function cellSubtitle() {
  const c = S.meta.counts;
  if (S.profileSel !== null) {
    const p = S.meta.profiles.find(p => p.id === S.profileSel);
    let s = `profile ${S.profileSel} — ${fmt(p ? p.n_cells : 0)} cells ` +
            `in ${fmt(p ? p.n_clone_units : 0)} clone units, ` +
            `against ${fmt(c.cells)} as context`;
    if (S.cellMode === 'celltype') s += ', coloured by cell type';
    return s;
  }
  let s = `${fmt(c.cells)} cells`;
  if (S.cellMode === 'profile') {
    const assigned = S.meta.profiles.reduce((a, p) => a + p.n_cells, 0);
    s += ` — ${fmt(assigned)} with a profile, ${fmt(c.cells - assigned)} without`;
  } else if (S.cellMode === 'celltype') {
    const ann = S.meta.cell_types.reduce((a, t) => a + t.n_cells, 0);
    s += ` — ${S.meta.cell_types.length} cell types`;
    if (ann < c.cells) s += `, ${fmt(c.cells - ann)} unannotated`;
  }
  if (S.isolate) s += ' — one class isolated';
  return s;
}

/* The MANAscore control exists only where the analysis has been run. It is hidden, not
 * disabled, on datasets without it -- a greyed-out control invites the question "why is
 * this off?", and the answer for Wang and Mathew is simply that no such data exists. */
function fillManaToggle() {
  const wrap = $('mana-wrap'), box = $('mana-toggle'), lab = $('mana-label');
  const ms = S.info.manascore;
  if (!ms) { wrap.hidden = true; box.checked = false; return; }
  wrap.hidden = false;
  box.checked = S.mana;
  lab.textContent = 'MANAscore';
  wrap.title =
    `${fmt(ms.cells_high)} of ${fmt(ms.cells_scored)} scored cells are MANAscore-high; ` +
    `${fmt(ms.clone_units_ptrc)} of ${fmt(ms.clone_units_called)} clone units are pTRC` +
    (ms.clone_units_no_cutoff
      ? `. ${fmt(ms.clone_units_no_cutoff)} clone units belong to a patient with no ` +
        `usable cutoff and could not be called either way — those get a dashed grey ring.`
      : '') +
    ` MANAscore is independent of Clonotrace, so agreement between the ring and the ` +
    `profile is evidence rather than restatement. Cells outside the scored set were ` +
    `never scored; they are not negatives.`;
}

/* The clone panel's layout choice.
 *
 * Hidden unless the dataset actually has more than one layout, so four of the five
 * datasets look exactly as they did. Switching re-fetches the clone positions and, if a
 * clone or a patient selection is open, that selection's frames too -- the arrows and
 * per-visit positions are coordinates and have to move with the layout. What is NOT
 * refetched is anything about identity: profiles, counts and the selection itself are
 * the same objects either way, which is why the current view survives the switch instead
 * of resetting to "all clones".
 */
function fillLayoutPick() {
  const sel = $('clone-layout');
  const wrap = $('clone-layout-wrap');
  const L = (S.meta && S.meta.clone_layouts) || [];
  sel.innerHTML = '';
  if (L.length < 2) { wrap.hidden = true; return; }
  for (const l of L) sel.appendChild(el('option', null, l.label)).value = l.id;
  sel.value = S.layout;
  wrap.hidden = false;
}

async function setLayout(id) {
  if (id === S.layout) return;
  const prev = S.layout;
  S.layout = id;
  try {
    S.clones = await api(`/api/clones/${S.dataset}?layout=${id}`);
    if (S.clone) {
      S.clone = await api(`/api/lineage/${S.dataset}/${S.clone.lineage}?layout=${id}`);
    } else if (S.set && S.patientsSel && S.patientsSel.length) {
      const q = S.patientsSel.map(p => `patient=${p}`).join('&');
      S.set = await api(`/api/set/${S.dataset}?${q}&layout=${id}`);
    }
  } catch (e) {
    S.layout = prev;
    $('clone-layout').value = prev;
    toast('could not switch layout: ' + e.message);
    return;
  }
  refresh();
}

function fillProfilePick() {
  const sel = $('profile-pick');
  sel.innerHTML = '';
  const all = document.createElement('option');
  all.value = ''; all.textContent = `all (${S.meta.profiles.length})`;
  sel.appendChild(all);
  for (const p of S.meta.profiles) {
    const o = document.createElement('option');
    o.value = String(p.id);
    o.textContent = `${p.id} — ${fmt(p.n_clone_units)} clone units`;
    sel.appendChild(o);
  }
  sel.value = S.profileSel === null ? '' : String(S.profileSel);
}

/* ---- legend, which is also the isolate control --------------------------- */
function renderLegend() {
  const wrap = $('legend'), title = $('legend-title'), qual = $('legend-quality');
  wrap.innerHTML = '';
  const byProfile = S.cellMode === 'profile' || S.cellMode === 'plain';
  const kind = byProfile ? 'profile' : 'celltype';
  title.textContent = byProfile ? 'Clone profiles' : 'Cell types';

  const q = byProfile ? S.info.profile_palette_quality : S.info.cell_type_palette_quality;
  qual.textContent = paletteNote(q, byProfile ? S.meta.profiles.length
                                              : S.meta.cell_types.length);

  const items = byProfile
    ? S.meta.profiles.map(p => ({ key: p.id, label: `profile ${p.id}`,
        colour: S.info.profile_colours[p.id], n: p.n_clone_units, unit: 'clone units' }))
    : S.meta.cell_types.map((t, i) => ({ key: i, label: t.name,
        colour: S.info.cell_type_colours[i], n: t.n_cells, unit: 'cells' }));

  for (const it of items) {
    const b = document.createElement('button');
    b.className = 'key';
    const on = S.isolate && S.isolate.kind === kind && S.isolate.value === it.key;
    b.setAttribute('aria-pressed', on ? 'true' : 'false');
    if (S.isolate && !on) b.classList.add('dim');
    b.innerHTML = `<span class="swatch" style="background:${it.colour}"></span>` +
                  `<span>${it.label}</span><span class="n">${fmt(it.n)}</span>`;
    b.title = `${it.label} — ${fmt(it.n)} ${it.unit}. Click to show only this one.`;
    b.onclick = () => {
      S.isolate = on ? null : { kind, value: it.key };
      refresh();
    };
    wrap.appendChild(b);
  }
  $('legend-clear').hidden = !S.isolate;
  if (S.mana && S.info.manascore) {
    const ms = S.info.manascore;
    const k = document.createElement('span');
    k.className = 'key mana-key';
    k.innerHTML =
      `<span class="swatch" style="background:${S.info.colours.mana}"></span>` +
      `<span>MANAscore-high cell / pTRC clone</span>` +
      `<span class="n">${fmt(ms.cells_high)} cells · ${fmt(ms.clone_units_ptrc)} clones</span>`;
    k.title = $('mana-wrap').title;
    wrap.appendChild(k);
  }
}

function paletteNote(q, n) {
  if (!q) return '';
  const cvd = q.worst_cvd, nm = q.worst_normal;
  const bits = [`${n} classes; worst pair ΔE ${cvd.toFixed(1)} under simulated ` +
                `colour-vision deficiency, ${nm.toFixed(1)} under normal vision`];
  if (q.meets_normal_floor === false || nm < 15) {
    bits.push('below the ΔE 15 floor — colour alone cannot separate these; ' +
              'use the labels, or click a class to isolate it');
  }
  return bits.join('. ') + '.';
}

/* ---- hover --------------------------------------------------------------- */
function wireHover() {
  const plot = $('clone-plot'), tip = $('clone-tip');
  plot.onmousemove = (e) => {
    if (!S.clones || !S.cloneFit) return;
    const r = plot.getBoundingClientRect();
    const px = e.clientX - r.left, py = e.clientY - r.top;
    const f = S.cloneFit, C = S.clones;
    let best = -1, bd = 12 * 12;
    for (let i = 0; i < C.x.length; i++) {
      const dx = f.X(C.x[i]) - px, dy = f.Y(C.y[i]) - py;
      const d = dx * dx + dy * dy;
      if (d < bd) { bd = d; best = i; }
    }
    if (best < 0) { tip.hidden = true; return; }
    const pid = S.info.patients[C.patient[best]];
    const v = S.meta.visits[C.visit[best]];
    tip.textContent =
      `${pid} · ${v.label}\nprofile ${C.profile[best]} · ${fmt(C.n_cells[best])} cells`;
    tip.style.left = px + 'px'; tip.style.top = py + 'px';
    tip.hidden = false;
  };
  plot.onmouseleave = () => { $('clone-tip').hidden = true; };
}

/* ---- orchestration ------------------------------------------------------- */
let centToken = 0;
async function refreshCentroids() {
  const my = ++centToken;
  const q = S.profileSel === null ? '' : `?profile=${S.profileSel}`;
  try {
    const c = await api(`/api/centroids/${S.dataset}${q}`);
    if (my !== centToken) return;
    S.centroids = c;
  } catch (e) {
    S.centroids = S.info.centroids;
  }
  drawCellLabels();
}

function refresh() {
  drawClones();
  updateCellBg();
  renderLegend();
  refreshCentroids();
  renderTimebar();
  renderCaptions();
  const a = active();
  if (a) {
    const f = currentFrame();
    if (f) { $('clone-sub').textContent = f.label; $('cell-sub').textContent = f.label; }
    if (a.kind === 'clone') ensureLineageRows().then(renderTimebar);
  }
  broadcastState();
}

async function openDataset(name) {
  S.dataset = name;
  $('dataset-line').textContent = 'loading…';
  S.layout = 'default';
  const [info, clones] = await Promise.all([
    api(`/api/dataset/${name}`), api(`/api/clones/${name}?layout=default`),
  ]);
  S.info = info; S.meta = info.meta; S.clones = clones;
  S.isolate = null; S.profileSel = null; S.centroids = info.centroids;
  S.patientSel = null; S.clone = null; S.frame = 0; S.lineageRows = null;
  S.mana = false;
  stopPlay();
  S.axis = 'time'; S.fixed = -1;
  S.patientsSel = []; S.lineagesSel = [];
  fillManaToggle();
  fillLayoutPick();
  fillAxisPick();
  fillProfilePick();
  const c = S.meta.counts;
  $('dataset-line').textContent =
    `${S.meta.label} — ${fmt(c.clone_units)} clone units over ${fmt(c.lineages)} ` +
    `lineages · ${fmt(c.patients)} patients · ${fmt(c.cells)} cells · ` +
    `${S.meta.profiles.length} profiles`;
  $('dataset-line').title = S.meta.notes || '';
  $('clone-note').textContent =
    `clone unit = one lineage in one patient at one visit` +
    (S.meta.expanded_min ? `, ≥${S.meta.expanded_min} cells` : '');
  /* Keeping the chosen dataset in the address bar makes the served viewer linkable. It
   * cannot be done offline: a file:// document has an opaque origin, and replaceState
   * refuses a URL that is not same-origin with the document -- which no URL is, when the
   * origin is "null". It throws rather than being ignored, and an uncaught throw here
   * would abort openDataset() before refresh(), leaving the page loaded and blank. The
   * address bar is a convenience; the guard is not. */
  try {
    const url = new URL(location);
    url.searchParams.set('dataset', name);
    history.replaceState({}, '', url);
  } catch (e) { /* file:// — the URL simply does not track the dataset */ }
  refresh();
}

async function boot() {
  const { datasets } = await api('/api/datasets');
  const pick = $('dataset-pick');
  if (!datasets.length) {
    $('dataset-line').textContent =
      'No prepared datasets. Run: python -m clonoview.prepare <name>';
    return;
  }
  for (const d of datasets) {
    const o = document.createElement('option');
    o.value = d.name; o.textContent = d.label;
    pick.appendChild(o);
  }
  pick.onchange = () => openDataset(pick.value).catch(e => toast(e.message));

  const want = new URL(location).searchParams.get('dataset');
  const start = datasets.find(d => d.name === want) ? want : datasets[0].name;
  pick.value = start;
  await openDataset(start);

  $('clone-colour').onchange = (e) => { S.cloneColour = e.target.checked; drawClones(); };
  $('clone-labels').onchange = (e) => { S.cloneLabels = e.target.checked; drawClones(); };
  $('clone-layout').onchange = (e) => { setLayout(e.target.value); };
  $('tb-axis').onchange = (e) => {
    // Switching axis resets what is held: the old value indexes the other axis's steps
    // and would silently mean something different. -1 ("all positions") is valid only on
    // the time axis; the tissue axis needs a real visit, and the server picks the busiest
    // one when given -1.
    setAxis(e.target.value, -1);
  };
  $('tb-fixed').onchange = (e) => { setAxis(S.axis, Number(e.target.value)); };
  $('cell-mode').onchange = (e) => { S.cellMode = e.target.value; S.isolate = null; refresh(); };
  $('mana-toggle').onchange = (e) => { S.mana = e.target.checked; refresh(); };
  $('profile-pick').onchange = (e) => {
    S.profileSel = e.target.value === '' ? null : Number(e.target.value);
    S.isolate = null;
    refresh();
  };
  $('cell-labels').onchange = (e) => { S.cellLabels = e.target.checked; drawCellLabels(); };
  $('legend-clear').onclick = () => { S.isolate = null; refresh(); };
  // The tables open in a SECOND TAB so they can sit beside the viewer; selections
  // travel back over a BroadcastChannel. See openTables().
  $('btn-patients').onclick = () => openTables('patients');
  $('btn-clones').onclick = () => openTables('clones');
  $('tb-close').onclick = clearSelection;
  $('tb-prev').onclick = () => { stopPlay(); if (S.frame > 0) { S.frame--; drawFrame(); } };
  $('tb-next').onclick = () => {
    stopPlay();
    if (S.frame < S.clone.frames.length - 1) { S.frame++; drawFrame(); }
  };
  $('tb-play').onclick = togglePlay;
  $('tb-trail').onchange = drawFrame;
  $('tb-contour').onchange = drawFrame;
  $('tb-prev-clone').onclick = () => stepClone(-1).catch(e => toast(e.message));
  $('tb-next-clone').onclick = () => stepClone(1).catch(e => toast(e.message));

  // Keyboard: arrows step time, space plays, Escape backs out one level. Cheap to add
  // and the difference between clicking through a time course and actually watching it.
  window.addEventListener('keydown', (e) => {
    if (e.target.matches('input, select, textarea')) return;
    if (!active()) return;
    if (e.key === 'ArrowLeft') { stopPlay(); if (S.frame > 0) { S.frame--; drawFrame(); } }
    else if (e.key === 'ArrowRight') {
      stopPlay();
      if (S.frame < S.clone.frames.length - 1) { S.frame++; drawFrame(); }
    } else if (e.key === ' ') { e.preventDefault(); togglePlay(); }
    else if (e.key === 'Escape') clearSelection();
    else return;
  });

  // Clicking a clone point opens that clone -- one click, as asked, same as the table.
  $('clone-plot').onclick = (e) => {
    if (!S.clones || !S.cloneFit) return;
    const r = $('clone-plot').getBoundingClientRect();
    const px = e.clientX - r.left, py = e.clientY - r.top;
    const f = S.cloneFit, C = S.clones;
    let best = -1, bd = 12 * 12;
    for (let i = 0; i < C.x.length; i++) {
      const dx = f.X(C.x[i]) - px, dy = f.Y(C.y[i]) - py;
      const d = dx * dx + dy * dy;
      if (d < bd) { bd = d; best = i; }
    }
    if (best < 0) return;
    // Shift (or ctrl/cmd) accumulates, matching the patient table's gesture, so a
    // handful of clones can be built up by clicking them on the map.
    if (e.shiftKey || e.ctrlKey || e.metaKey) toggleLineage(C.lineage[best]);
    else selectClone(C.lineage[best]);
  };

  wireHover();
  let t = null;
  window.addEventListener('resize', () => {
    clearTimeout(t);
    t = setTimeout(() => { drawClones(); updateCellBg(); }, 120);
  });
}

boot().catch(e => { $('dataset-line').textContent = 'error: ' + e.message; });

/* =========================================================================
 * PATIENT TABLE, CLONE TABLE, AND THE PER-CLONE VIEW
 * ========================================================================= */

/* =========================================================================
 * ONE SELECTION MODEL FOR BOTH VIEWS
 *
 * A single clone and a set of clones are the same object as far as drawing is
 * concerned: an ordered list of frames, one per visit, each carrying positions on the
 * clone map, arrows from the previous visit, cells, and a density. Keeping them on one
 * code path is what stops the set view from quietly diverging from the clone view --
 * which is exactly how the cell panel ended up ignoring a patient selection.
 * ========================================================================= */
const SVGNS = 'http://www.w3.org/2000/svg';
function mk(tag, attrs) {
  const e = document.createElementNS(SVGNS, tag);
  for (const k in attrs) e.setAttribute(k, attrs[k]);
  return e;
}

function active() {
  if (S.clone) {
    return { kind: 'clone', frames: S.clone.frames, label: S.clone.label,
             nCells: S.clone.n_cells };
  }
  if (S.set) {
    const who = S.set.patient_ids.length ? S.set.patient_ids.join(', ') : 'selection';
    return { kind: 'set', frames: S.set.frames, label: who, nCells: S.set.n_cells };
  }
  return null;
}

function currentFrame() {
  const a = active();
  if (!a) return null;
  return a.frames[Math.min(S.frame, a.frames.length - 1)] || null;
}

/* ---- LEFT PANEL: positions and arrows ------------------------------------ */
function drawTrack(g, f) {
  const a = active();
  if (a.kind === 'clone') return drawCloneTrack(g, f);
  return drawSetTrack(g, f);
}

/* One clone: its position at each visit, joined in time order. Frames with no position
 * are skipped by the arrows but not hidden -- the right panel still draws their cells,
 * and an arrow that jumps over one is dashed, because it spans a gap not a step. */
function drawCloneTrack(g, f) {
  const fr = S.clone.frames;
  const withPos = fr.map((x, i) => [x, i]).filter(([x]) => x.xy);
  const trail = $('tb-trail').checked;

  for (let k = 1; k < withPos.length; k++) {
    const [aF, ai] = withPos[k - 1], [bF, bi] = withPos[k];
    if (!trail && bi > S.frame) break;
    arrow(g, f.X(aF.xy[0]), f.Y(aF.xy[1]), f.X(bF.xy[0]), f.Y(bF.xy[1]),
          bi === S.frame, bi - ai > 1, trail);
  }
  for (const [x, i] of withPos) {
    if (!trail && i > S.frame) continue;
    const cur = i === S.frame;
    dot(g, f.X(x.xy[0]), f.Y(x.xy[1]), cur ? 7 : 4.2,
        S.info.profile_colours[x.profile] || '#6d6b64', cur);
  }
}

/* A set of clones: every clone unit present at the current visit, and every arrow
 * arriving at it. With "Trail" on, earlier visits stay faint so the accumulated
 * movement is visible; with it off you see one step at a time. */
function drawSetTrack(g, f) {
  const trail = $('tb-trail').checked;
  const frames = S.set.frames;
  for (let i = 0; i <= S.frame; i++) {
    const fr = frames[i];
    const cur = i === S.frame;
    if (!cur && !trail) continue;
    for (const ar of fr.arrows) {
      arrow(g, f.X(ar.x1), f.Y(ar.y1), f.X(ar.x2), f.Y(ar.y2), cur, ar.skipped, trail);
    }
  }
  for (let i = 0; i <= S.frame; i++) {
    const fr = frames[i];
    const cur = i === S.frame;
    if (!cur && !trail) continue;
    const u = fr.units;
    for (let k = 0; k < u.x.length; k++) {
      dot(g, f.X(u.x[k]), f.Y(u.y[k]), cur ? 4.6 : 3.0,
          S.info.profile_colours[u.profile[k]] || '#6d6b64', cur,
          markerFor(u.lineage[k]));
    }
  }
}

function arrow(g, x1, y1, x2, y2, active_, skipped, trail) {
  g.save();
  g.strokeStyle = active_ ? '#0b0b0b' : '#8a8880';
  g.lineWidth = active_ ? 1.8 : 1.0;
  g.globalAlpha = active_ ? 0.95 : (trail ? 0.4 : 0.95);
  if (skipped) g.setLineDash([4, 3]);
  g.beginPath(); g.moveTo(x1, y1); g.lineTo(x2, y2); g.stroke();
  const ang = Math.atan2(y2 - y1, x2 - x1), L = active_ ? 8 : 6;
  g.setLineDash([]);
  g.beginPath();
  g.moveTo(x2, y2);
  g.lineTo(x2 - L * Math.cos(ang - 0.4), y2 - L * Math.sin(ang - 0.4));
  g.lineTo(x2 - L * Math.cos(ang + 0.4), y2 - L * Math.sin(ang + 0.4));
  g.closePath();
  g.fillStyle = g.strokeStyle; g.fill();
  g.restore();
}

/* ---- markers -------------------------------------------------------------
 * Shape is the second categorical channel: colour already means profile. The geometry
 * and the honest cap on how many clones can be told apart live in static/markers.js,
 * shared with the clone table so the two cannot disagree.
 */

/* Which shape a lineage gets, or null for "not individually marked". Driven by
 * `lineages_selected`, which the server sorts, so the mapping is stable across
 * reloads, across frames, and across the two panels. */
function markerFor(lineage) {
  const sel = (S.set && S.set.lineages_selected) || [];
  return MARKERS.shapeAt(sel.indexOf(lineage), sel.length);
}

function dot(g, px, py, r, fill, emphasised, shape) {
  MARKERS.path(g, shape || 'circle', px, py, r);
  g.fillStyle = fill;
  g.globalAlpha = emphasised ? 1 : 0.7;
  g.fill();
  g.globalAlpha = 1;
  g.lineWidth = emphasised ? 2.0 : 1.0;
  g.strokeStyle = emphasised ? '#0b0b0b' : '#fcfcfb';
  g.stroke();
}

/* ---- RIGHT PANEL: cells on canvas, contours on SVG ----------------------- */
function drawCellSelection() {
  const cv = $('cell-canvas'), plot = $('cell-plot');
  const w = plot.clientWidth, h = plot.clientHeight;
  const dpr = window.devicePixelRatio || 1;
  cv.width = Math.round(w * dpr); cv.height = Math.round(h * dpr);
  const g = cv.getContext('2d');
  g.setTransform(dpr, 0, 0, dpr, 0, 0);
  g.clearRect(0, 0, w, h);

  const a = active();
  const fr = currentFrame();
  if (!a || !fr) return;
  const ext = S.meta.extent.cell;
  const fit = fitBox(ext, w, h);
  const colour = a.kind === 'clone'
    ? (fr.below_threshold ? '#0b0b0b'
                          : (S.info.profile_colours[fr.profile] || '#0b0b0b'))
    : '#0b0b0b';

  const xs = fr.cells.x, ys = fr.cells.y;
  const r = a.kind === 'clone' ? 2.6 : 1.5;
  g.fillStyle = colour;
  g.strokeStyle = '#fcfcfb';
  g.lineWidth = a.kind === 'clone' ? 1.1 : 0;
  for (let i = 0; i < xs.length; i++) {
    g.beginPath(); g.arc(fit.X(xs[i]), fit.Y(ys[i]), r, 0, 7); g.fill();
    if (a.kind === 'clone') g.stroke();
  }

  // Contours go on the SVG, in data units, so they scale with the viewBox and sit
  // above the cells. Cells are the observation; the ring is a model of them, so the
  // ring is translucent and never hides a cell.
  if (!$('tb-contour').checked) return;
  const d = fr.density;
  if (!d || !d.rings) return;
  const svg = $('cell-overlay');
  const span = ext[1] - ext[0];
  d.rings.forEach((polys, li) => {
    const frac = d.levels[li];
    for (const flat of polys) {
      let path = '';
      for (let i = 0; i < flat.length; i += 2) {
        path += (i ? 'L' : 'M') + flat[i] + ' ' + (-flat[i + 1]) + ' ';
      }
      svg.appendChild(mk('path', {
        d: path + 'Z', fill: colour, 'fill-opacity': 0.07,
        stroke: colour, 'stroke-opacity': 0.35 + 0.2 * li,
        'stroke-width': span / (frac > 0.7 ? 900 : 650),
        'stroke-linejoin': 'round',
      }));
    }
  });
}

function drawFrame() {
  drawClones();
  drawCellLabels();          // clears the SVG, then re-adds labels and contours
  renderTimebar();
  renderCaptions();
  const f = currentFrame();
  if (f) { $('clone-sub').textContent = f.label; $('cell-sub').textContent = f.label; }
}

/* ---- time bar, shared by both views -------------------------------------- */
function renderTimebar() {
  const bar = $('timebar');
  const a = active();
  if (!a) { bar.hidden = true; return; }
  bar.hidden = false;

  const unit = S.axis === 'tissue' ? 'anatomic position' : 'timepoint';
  $('tb-what').textContent = a.kind === 'clone'
    ? `${a.label}  ·  ${fmt(a.nCells)} cells across ${a.frames.length} ${unit}` +
      `${a.frames.length === 1 ? '' : 's'}`
    : `${a.label}  ·  ${fmt(S.set.n_lineages)} clones, ${fmt(S.set.n_units)} clone ` +
      `units, ${fmt(a.nCells)} cells`;
  fillAxisPick();

  /* On the anatomic axis the clone map deliberately does not move: Clonotrace placed
   * clone units per visit, so there is no per-site position to step to. Saying so is the
   * difference between a reader seeing a limitation and a reader seeing a bug. */
  const held = S.axis === 'tissue'
    ? (S.meta.visits[S.fixed] || {}).label
    : (S.fixed >= 0 ? (S.meta.tissues || [])[S.fixed] &&
       (S.meta.tissues || [])[S.fixed].label : null);
  $('clone-note').textContent = S.axis === 'tissue'
    ? `stepping through anatomic positions at ${held} — clone units are placed per ` +
      `visit, so the clone map holds still and only the cells move`
    : (held ? `clone unit = one lineage in one patient at one visit` +
              (S.meta.expanded_min ? `, ≥${S.meta.expanded_min} cells` : '') +
              ` · cells restricted to ${held}`
            : `clone unit = one lineage in one patient at one visit` +
              (S.meta.expanded_min ? `, ≥${S.meta.expanded_min} cells` : ''));

  const steps = $('tb-steps');
  steps.innerHTML = '';
  a.frames.forEach((fr, i) => {
    const noPos = a.kind === 'clone' ? fr.below_threshold : fr.n_units === 0;
    const b = el('button', 'step' + (noPos ? ' nopos' : ''), fr.code);
    b.setAttribute('aria-pressed', i === S.frame ? 'true' : 'false');
    if (a.kind === 'clone' && !noPos) {
      b.style.borderColor = S.info.profile_colours[fr.profile];
      if (i === S.frame) b.style.background = S.info.profile_colours[fr.profile];
    }
    b.title = a.kind === 'clone'
      ? `${fr.label} — ${fmt(fr.n_cells)} cells` +
        (noPos ? `, below the ${S.meta.expanded_min}-cell threshold, so no ` +
                 `clone-map position` : `, profile ${fr.profile}`)
      : `${fr.label} — ${fmt(fr.n_units)} clone units, ${fmt(fr.n_cells)} cells`;
    b.onclick = () => { stopPlay(); S.frame = i; drawFrame(); };
    steps.appendChild(b);
  });
  $('tb-prev').disabled = S.frame <= 0;
  $('tb-next').disabled = S.frame >= a.frames.length - 1;

  const single = a.kind === 'clone';
  const idx = single && S.lineageRows
    ? S.lineageRows.findIndex(r => r.lineage === S.clone.lineage) : -1;
  $('tb-prev-clone').hidden = !single;
  $('tb-next-clone').hidden = !single;
  $('tb-prev-clone').disabled = idx <= 0;
  $('tb-next-clone').disabled = idx < 0 || idx >= (S.lineageRows || []).length - 1;
}

/* The shape key for a multi-clone selection. Deliberately says nothing about colour:
 * colour is the profile, in this panel as in every other, and a selection of three
 * clones does not change that. */
function renderMarkerKey() {
  const box = $('clone-marker-key');
  if (!box) return;
  const sel = (S.set && S.set.lineages_selected) || [];
  const names = (S.set && S.set.lineages_selected_labels) || [];
  const name = (i) => names[i] || ('clone ' + sel[i]);
  if (sel.length < 2) { box.hidden = true; box.innerHTML = ''; return; }
  box.hidden = false;
  if (sel.length > MARKERS.MAX) {
    box.innerHTML = `<span class="note">${fmt(sel.length)} clones selected — too many ` +
      `to mark apart, so all are drawn as circles. Select ${MARKERS.MAX} or fewer to ` +
      `give each its own shape.</span>`;
    return;
  }
  box.innerHTML = '<span class="note">Shape = clone · colour = profile</span>' +
    sel.map((l, i) => `<span class="mk-item" title="${name(i)}">` +
      `${MARKERS.glyph(MARKERS.SHAPES[i])}<span>${name(i)}</span></span>`).join('');
}

function renderCaptions() {
  renderMarkerKey();
  const cc = $('clone-caption'), xc = $('cell-caption');
  const a = active(), f = currentFrame();
  if (!a || !f) { cc.hidden = true; xc.hidden = true; return; }
  cc.hidden = false; xc.hidden = false;
  cc.innerHTML = ''; xc.innerHTML = '';

  cc.appendChild(el('span', 'id', a.label));
  if (a.kind === 'clone') {
    cc.appendChild(el('span', 'n', `${fmt(f.n_cells)} cells`));
    cc.appendChild(f.below_threshold
      ? el('span', 'warn', `no position — under ${S.meta.expanded_min} cells here`)
      : el('span', 'n', `profile ${f.profile}`));
  } else {
    cc.appendChild(el('span', 'n', `${fmt(f.n_units)} clone units here`));
    if (f.n_units === 0 && f.n_cells === 0) {
      cc.appendChild(el('span', 'warn', 'not sampled at this visit'));
    }
  }

  xc.appendChild(el('span', 'id', a.label));
  xc.appendChild(el('span', 'n', `${fmt(f.n_cells)} cells`));
  if (f.subsampled) {
    xc.appendChild(el('span', 'n',
      `showing ${fmt(f.cells_shown)} of them; contours use all ${fmt(f.n_cells)}`));
  }
  const d = f.density || {};
  if (d.reason) xc.appendChild(el('span', 'n', d.reason));
  else if (d.levels && d.levels.length) {
    xc.appendChild(el('span', 'n', 'rings enclose ' +
      d.levels.map(l => `${Math.round(l * 100)}%`).join(', ') + ' of these cells'));
  }
}

function el(tag, cls, text) {
  const e = document.createElement(tag);
  if (cls) e.className = cls;
  if (text !== undefined) e.textContent = text;
  return e;
}

/* ---- selecting ----------------------------------------------------------- */
async function selectClone(lineage) {
  stopPlay();
  try {
    S.clone = await api(`/api/lineage/${S.dataset}/${lineage}` +
                        `?layout=${S.layout}&axis=${S.axis}&fixed=${S.fixed}`);
  } catch (e) { toast('could not load that clone: ' + e.message); return; }
  if (S.clone.axis) { S.axis = S.clone.axis; S.fixed = S.clone.fixed; }
  S.set = null;
  S.frame = 0;
  S.patientSel = S.clone.patient;
  S.patientsSel = [];
  S.lineagesSel = [lineage];
  S.profileSel = null;
  fillProfilePick();
  fillAxisPick();
  refresh();
}

/* Any selection of more than one thing -- several patients, several clones, or a mix of
 * both. /api/set unions `patient=` and `lineage=`, so the two kinds live in one request
 * and one view; nothing here needs to know which kind the user clicked first. */
async function selectMany(patients, lineages) {
  stopPlay();
  patients = patients || []; lineages = lineages || [];
  if (!patients.length && !lineages.length) { clearSelection(); return; }
  if (!patients.length && lineages.length === 1) { return selectClone(lineages[0]); }
  const q = [
    ...patients.map(p => `patient=${p}`),
    ...lineages.map(l => `lineage=${l}`),
    `layout=${S.layout}`, `axis=${S.axis}`, `fixed=${S.fixed}`,
  ].join('&');
  try { S.set = await api(`/api/set/${S.dataset}?${q}`); }
  catch (e) { toast('could not load that selection: ' + e.message); return; }
  if (S.set.axis) { S.axis = S.set.axis; S.fixed = S.set.fixed; }
  S.clone = null;
  S.frame = 0;
  S.patientsSel = patients;
  S.lineagesSel = lineages;
  S.patientSel = patients.length === 1 ? patients[0] : null;
  S.lineageRows = null;
  S.profileSel = null;
  fillProfilePick();
  fillAxisPick();
  refresh();
}

/* Kept for the existing callers and the cross-tab message of the same name. */
async function selectPatients(patients) {
  return selectMany(patients, []);
}

/* Add or remove one clone from whatever is currently selected. This is what makes
 * shift-click on the clone map and in the clone table accumulate rather than replace. */
function toggleLineage(lineage) {
  const have = (S.lineagesSel || []).includes(lineage);
  const lins = have ? S.lineagesSel.filter(l => l !== lineage)
                    : [...(S.lineagesSel || []), lineage];
  return selectMany(S.patientsSel || [], lins);
}

function clearSelection() {
  stopPlay();
  S.clone = null; S.set = null; S.patientSel = null;
  S.patientsSel = []; S.lineagesSel = [];
  S.frame = 0;
  refresh();
}

/* ---- the two axes -------------------------------------------------------
 * Time and anatomic position are both per-cell covariates, so either can be the one
 * that steps. They are NOT symmetric on the clone map: Clonotrace placed clone units at
 * (patient, visit, receptor), so a clone has a position per visit and none per tissue.
 * Stepping through anatomic positions therefore holds the left panel still and varies
 * only the right, and the note under the panel says so rather than leaving the reader to
 * wonder why the points stopped moving.
 */
function fillAxisPick() {
  const axes = (S.info && S.info.axes) || [];
  const aw = $('tb-axis-wrap'), fw = $('tb-fixed-wrap');
  if (axes.length < 2) { aw.hidden = true; fw.hidden = true; return; }
  const sel = $('tb-axis');
  sel.innerHTML = '';
  for (const a of axes) sel.appendChild(el('option', null, a.label)).value = a.id;
  sel.value = S.axis;
  aw.hidden = false;

  // The "holding" list is the OTHER axis. On the time axis it may be "all positions";
  // on the tissue axis a visit must be named, because a clone has one clone-map position
  // per visit and stepping through sites has to hold one of them.
  const other = axes.find(a => a.id !== S.axis);
  const f = $('tb-fixed');
  f.innerHTML = '';
  if (S.axis === 'time') {
    f.appendChild(el('option', null, 'all positions')).value = '-1';
  }
  other.steps.forEach((st, i) => {
    const o = el('option', null, st.label + (st.n_cells ? ` (${fmt(st.n_cells)})` : ''));
    o.value = String(i);
    f.appendChild(o);
  });
  f.value = String(S.fixed);
  if (!f.value || f.selectedIndex < 0) f.value = f.options[0].value;
  fw.hidden = false;
}

async function setAxis(axis, fixed) {
  S.axis = axis;
  S.fixed = fixed;
  S.frame = 0;
  if (S.clone) await selectClone(S.clone.lineage);
  else if (S.set) await selectMany(S.patientsSel || [], S.lineagesSel || []);
  else fillAxisPick();
}

/* ---- playback ------------------------------------------------------------ */
let playTimer = null;
function stopPlay() {
  if (playTimer) { clearInterval(playTimer); playTimer = null; }
  const b = $('tb-play');
  if (b) { b.textContent = '▶ Play'; b.setAttribute('aria-pressed', 'false'); }
}
function togglePlay() {
  if (playTimer) { stopPlay(); return; }
  const a = active();
  if (!a || a.frames.length < 2) return;
  const b = $('tb-play');
  b.textContent = '❚❚ Pause'; b.setAttribute('aria-pressed', 'true');
  playTimer = setInterval(() => {
    S.frame = (S.frame + 1) % active().frames.length;
    drawFrame();
  }, 1100);
}

async function stepClone(delta) {
  if (!S.lineageRows || !S.clone) return;
  const i = S.lineageRows.findIndex(r => r.lineage === S.clone.lineage);
  const j = i + delta;
  if (j < 0 || j >= S.lineageRows.length) return;
  await selectClone(S.lineageRows[j].lineage);
}

async function ensureLineageRows() {
  if (S.lineageRows && S.clone
      && S.lineageRows.some(r => r.lineage === S.clone.lineage)) return;
  const q = S.patientSel === null ? '' : `?patient=${S.patientSel}`;
  try {
    const d = await api(`/api/lineages/${S.dataset}${q}`);
    S.lineageRows = d.rows;
  } catch (e) { /* stepping simply stays disabled */ }
}

/* ---- the tables, in a second tab -----------------------------------------
 * Opened with a fixed window name so repeated clicks reuse the same tab rather than
 * piling up. Selections come back over the LINK (see link.js); the viewer also publishes
 * its own state so the table can show which row is currently displayed.
 */
function openTables(which) {
  // Relative, not "/static/table.html": offline there is no server and no site root, and
  // an absolute path would resolve to the filesystem root and open nothing.
  const base = LOCAL() ? 'table.html' : '/static/table.html';
  const url = `${base}?dataset=${encodeURIComponent(S.dataset)}&view=${which}`;
  const w = window.open(url, 'clonoview-tables');
  if (!w) { toast('Your browser blocked the new tab — allow pop-ups for this page.'); return; }
  LINK.attachChild(w);               // the handle IS the transport when running offline
  w.focus();
  setTimeout(broadcastState, 300);   // the new tab needs the current state on arrival
}

function broadcastState() {
  if (!S.meta) return;
  LINK.send({
    type: 'state', dataset: S.dataset,
    lineage: S.clone ? S.clone.lineage : null,
    patients: S.patientsSel || [],
    lineages: S.lineagesSel || [],
    axis: S.axis, fixed: S.fixed,
  });
}

LINK.onMessage((m) => {
  if (m.type === 'select') selectMany(m.patients || [], m.lineages || []);
  else if (m.type === 'selectClone') selectClone(m.lineage);
  else if (m.type === 'selectPatients') selectPatients(m.patients);
  else if (m.type === 'clear') clearSelection();
  else if (m.type === 'hello') broadcastState();
});
