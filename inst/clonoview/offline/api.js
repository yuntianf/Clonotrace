/* The server's JSON routes, reimplemented in the browser.
 *
 * This is a port of the `Api` class in clonoview/server.py, method for method, and the
 * comments there are the explanation -- what a clone unit is and is not, why centroids
 * are recomputed under a filter instead of reused, why the whole of a clone's time course
 * comes back at once. Only what is different here is written down.
 *
 * WHY IT KEEPS THE URL SHAPE. Every route returns exactly the object the served version
 * returns, and dispatch is on the same "/api/lineage/mathew/12" strings. The alternative
 * -- a clean local API and a rewritten front end -- would have been tidier and much
 * worse: app.js and table.js are 1,200 lines whose rendering has only ever been tested by
 * being looked at in a browser, and there is no JavaScript runtime on the machine this
 * was built on to test them any other way. Keeping the seam at api() means the served
 * version and the offline export run the SAME front-end code, so the served version stays
 * usable as the reference to compare against, and a rendering bug found in one is a bug
 * in both rather than a difference between two forks.
 *
 * WHAT IS DELIBERATELY NOT COPIED: the rounding. server.py rounds coordinates to 4 or 2
 * decimals on the way out, which exists to cut the size of a JSON payload crossing an SSH
 * tunnel. There is no tunnel here and no payload -- these objects are handed straight to
 * the caller -- so the values are passed at full precision. selftest.html therefore
 * compares against Python with a tolerance of half the rounding step it dropped, which is
 * a number with a reason rather than a fudge factor.
 */

window.CV = window.CV || {};

(function (CV) {
  'use strict';

  const MAX_CELLS_PER_FRAME = 12000;   // server.py Api.MAX_CELLS_PER_FRAME

  function median(vals) {
    if (!vals.length) return NaN;
    const a = Float64Array.from(vals).sort();
    const m = a.length;
    // np.median averages the two middle values on an even count. Taking the upper of
    // the two instead moves a direct label by up to half the gap between two cells,
    // which is invisible in the middle of a cloud and not invisible at the edge of one.
    return (m % 2) ? a[(m - 1) / 2] : (a[m / 2 - 1] + a[m / 2]) / 2;
  }

  /* numpy's default_rng(0).shuffle, which render.py uses to make occlusion unbiased, is
   * PCG64 and is not reproducible in JavaScript without porting the generator. The
   * requirement is a FIXED permutation, not numpy's particular one -- what must not
   * happen is drawing class by class, or reshuffling on every redraw so the same overlap
   * shows a different class each time the panel is touched. A small deterministic PRNG
   * satisfies that; the self-check compares the count of cells drawn per class, which is
   * the property that carries the meaning, not their order. */
  function mulberry32(seed) {
    let a = seed >>> 0;
    return function () {
      a = (a + 0x6D2B79F5) >>> 0;
      let t = a;
      t = Math.imul(t ^ (t >>> 15), t | 1);
      t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
  }

  function shuffleInPlace(idx, seed) {
    const rnd = mulberry32(seed === undefined ? 0 : seed);
    for (let i = idx.length - 1; i > 0; i--) {
      const j = Math.floor(rnd() * (i + 1));
      const t = idx[i]; idx[i] = idx[j]; idx[j] = t;
    }
    return idx;
  }

  /* ---- the routes --------------------------------------------------------- */

  function datasets() {
    return { datasets: window.CV_INDEX.datasets };
  }

  /* Mirrors Bundle.axes(): time always, anatomic position only where a cohort was
   * sampled at more than one site. */
  function tissuesOf(b) { return (b.meta.tissues || []).map((t) => t.code); }
  function axesOf(b) {
    const out = [{ id: 'time', label: 'Time',
                   steps: b.meta.visits.map((v) => ({ code: v.code, label: v.label })) }];
    const T = b.meta.tissues || [];
    if (T.length > 1) {
      out.push({ id: 'tissue', label: 'Anatomic position',
                 steps: T.map((t) => ({ code: t.code, label: t.label, n_cells: t.n_cells })) });
    }
    return out;
  }

  function dataset(b) {
    const nprof = b.profile_ids.length;
    const ctc = {};
    for (let i = 0; i < b.cell_type_colours.length; i++) ctc[String(i)] = b.cell_type_colours[i];
    const pc = {};
    for (const k in b.profile_colours) pc[String(k)] = b.profile_colours[k];
    return {
      meta: b.meta,
      patients: b.names.patients,
      cell_types: b.cell_types,
      profile_colours: pc,
      cell_type_colours: ctc,
      cell_type_palette_quality: b.colours.cell_type_palette_quality,
      profile_palette_quality: b.colours.profile_palette_quality,
      colours: b.colours.ui,
      manascore: b.meta.manascore === undefined ? null : b.meta.manascore,
      axes: axesOf(b),
      centroids: {
        profile: centroidsOf(b, b.c_profile, b.profile_ids, null, false),
        cell_type: centroidsOf(b, b.c_type, range(b.cell_types.length), null, false),
      },
    };
  }

  function range(n) { const a = []; for (let i = 0; i < n; i++) a.push(i); return a; }

  /* Median position per class, optionally over one profile's cells only, and -- when
   * `withCount` -- how many cells that is. Both callers are in server.py: the dataset
   * route wants [x, y] and the centroids route wants [x, y, n]. Same computation, two
   * payload shapes, kept together so a change to one cannot miss the other. */
  function centroidsOf(b, values, labels, keepProfile, withCount) {
    const out = {};
    for (const lab of labels) {
      const xs = [], ys = [];
      for (let i = 0; i < values.length; i++) {
        if (values[i] !== lab) continue;
        if (keepProfile !== null && b.c_profile[i] !== keepProfile) continue;
        xs.push(b.cxy[2 * i]); ys.push(b.cxy[2 * i + 1]);
      }
      if (!xs.length) continue;
      out[lab] = withCount ? [median(xs), median(ys), xs.length]
                           : [median(xs), median(ys)];
    }
    return out;
  }

  function centroids(b, q) {
    const p = q.get('profile');
    const prof = (p === null || p === '' || p === 'null') ? null : parseInt(p, 10);
    return {
      profile: centroidsOf(b, b.c_profile, b.profile_ids, prof, true),
      cell_type: centroidsOf(b, b.c_type, range(b.cell_types.length), prof, true),
      filtered_to_profile: prof,
    };
  }

  /* Clone coordinates under one layout. Mirrors Bundle.uxy_for on the Python side,
   * including its fallback: an unknown id gives the default rather than throwing, so a
   * stale link shows the map instead of an error. */
  function uxyOf(b, layout) {
    const id = layout || 'default';
    return (id !== 'default' && b['uxy_' + id]) ? b['uxy_' + id] : b.uxy;
  }
  function cloneExtentOf(b, layout) {
    const L = (b.meta.clone_layouts || []).filter((l) => l.id === (layout || 'default'));
    return L.length ? L[0].extent : b.meta.extent.clone;
  }

  function clones(b, q) {
    const n = b.n_units;
    const x = new Array(n), y = new Array(n);
    const layout = q ? (q.get('layout') || 'default') : 'default';
    const uxy = uxyOf(b, layout);
    for (let i = 0; i < n; i++) { x[i] = uxy[2 * i]; y[i] = uxy[2 * i + 1]; }
    return {
      x: x, y: y,
      profile: Array.from(b.u_profile),
      patient: Array.from(b.u_patient),
      visit: Array.from(b.u_visit),
      lineage: Array.from(b.u_lineage),
      n_cells: Array.from(b.u_n),
      ptrc: Array.from(b.u_ptrc),
      mana_usable: Array.from(b.u_mana_usable),
      layout: layout,
      extent: cloneExtentOf(b, layout),
    };
  }

  function lineages(b, q) {
    const want = q.getAll('patient').map((v) => parseInt(v, 10));
    const pats = want.length ? new Set(want) : null;
    // One pass over the clone units to group them by lineage, rather than scanning the
    // whole array once per lineage. Wang has 703 lineages and 1,514 units; the quadratic
    // form is what made the served version's table feel slow to open on it.
    const unitsBy = new Map();
    for (let k = 0; k < b.n_units; k++) {
      const li = b.u_lineage[k];
      if (!unitsBy.has(li)) unitsBy.set(li, []);
      unitsBy.get(li).push(k);
    }
    const rows = [];
    const L = b.names.lineages;
    for (let i = 0; i < L.length; i++) {
      const l = L[i];
      if (pats !== null && !pats.has(l.patient)) continue;
      const units = (unitsBy.get(i) || []).slice();
      units.sort((p, r) => b.u_visit[p] - b.u_visit[r]);
      const span = b.cellsOfLineage(i);
      const seen = new Set();
      for (let c = span[0]; c < span[1]; c++) {
        const v = b.c_visit[c];
        if (v >= 0) seen.add(v);
      }
      rows.push({
        lineage: i, label: l.label, patient: l.patient,
        receptor: l.receptor, n_cells: l.n_cells,
        n_units: units.length,
        visits: units.map((k) => b.u_visit[k]).sort((p, r) => p - r),
        profiles: units.map((k) => b.u_profile[k]),
        visits_with_cells: Array.from(seen).sort((p, r) => p - r),
      });
    }
    rows.sort((p, r) => (r.n_units - p.n_units) || (r.n_cells - p.n_cells));
    return { rows: rows, n: rows.length };
  }

  /* One frame per visit at which the lineage has any cell, plus any visit at which it was
   * a clone unit. bundle.py lineage_frames() is the original; the cell slice is turned
   * into the cells themselves here, as server.py's lineage route does. */
  function lineageFrames(b, li, layout, axis, fixed) {
    const uxy = uxyOf(b, layout);
    const span = b.cellsOfLineage(li);
    const a = span[0], z = span[1];
    const tissues = tissuesOf(b);
    const units = new Map();
    for (let k = 0; k < b.n_units; k++) {
      if (b.u_lineage[k] === li) units.set(b.u_visit[k], k);
    }
    /* Mirrors Bundle.lineage_frames, including the asymmetry: clone units exist per
     * visit and never per tissue, so the tissue axis carries the SAME unit on every
     * frame and only the cells move. */
    function mkFrame(vi, ti, k) {
      let first = -1, last = -1, count = 0;
      for (let c = a; c < z; c++) {
        if (b.c_visit[c] !== vi) continue;
        if (ti >= 0 && b.c_tissue[c] !== ti) continue;
        if (first < 0) first = c;
        last = c; count++;
      }
      const f = {
        visit: vi,
        code: axis === 'tissue' ? tissues[ti] : b.visits[vi],
        label: axis === 'tissue' ? tissues[ti] : b.meta.visits[vi].label,
        n_cells: count,
        cell_slice: count ? [first, last + 1] : null,
        clone_unit: k,
        xy: k !== null ? [uxy[2 * k], uxy[2 * k + 1]] : null,
        profile: k !== null ? b.u_profile[k] : null,
        below_threshold: k === null,
      };
      if (ti >= 0) { f.tissue = ti; f.tissue_label = tissues[ti]; }
      return f;
    }
    const frames = [];
    if (axis === 'tissue') {
      const vi = fixed;
      const k = units.has(vi) ? units.get(vi) : null;
      for (let ti = 0; ti < tissues.length; ti++) {
        const f = mkFrame(vi, ti, k);
        if (!f.n_cells) continue;
        f.label = tissues[ti]; f.code = tissues[ti];
        frames.push(f);
      }
      return frames;
    }
    const ti = (fixed >= 0 && fixed < tissues.length) ? fixed : -1;
    for (let vi = 0; vi < b.visits.length; vi++) {
      const k = units.has(vi) ? units.get(vi) : null;
      const f = mkFrame(vi, ti, k);
      if (!f.n_cells && k === null) continue;
      frames.push(f);
    }
    return frames;
  }

  function busiestVisit(b, pred) {
    const n = new Array(b.visits.length).fill(0);
    for (let c = 0; c < b.n_cells; c++) if (pred(c) && b.c_visit[c] >= 0) n[b.c_visit[c]]++;
    let best = 0;
    for (let i = 1; i < n.length; i++) if (n[i] > n[best]) best = i;
    return best;
  }

  function lineage(b, idx, q) {
    const i = parseInt(idx, 10);
    const l = b.names.lineages[i];
    const layout = q ? (q.get('layout') || 'default') : 'default';
    const axis = q ? (q.get('axis') || 'time') : 'time';
    let fixed = q ? parseInt(q.get('fixed') || '-1', 10) : -1;
    if (axis === 'tissue' && fixed < 0) {
      const span = b.cellsOfLineage(i);
      fixed = busiestVisit(b, (c) => c >= span[0] && c < span[1]);
    }
    const frames = lineageFrames(b, i, layout, axis, fixed);
    for (const f of frames) {
      if (f.cell_slice) {
        const [s, e] = f.cell_slice;
        const m = e - s;
        const x = new Array(m), y = new Array(m), ct = new Array(m), pr = new Array(m);
        for (let c = 0; c < m; c++) {
          x[c] = b.cxy[2 * (s + c)]; y[c] = b.cxy[2 * (s + c) + 1];
          ct[c] = b.c_type[s + c]; pr[c] = b.c_profile[s + c];
        }
        f.cells = { x: x, y: y, cell_type: ct, profile: pr };
      } else {
        f.cells = { x: [], y: [], cell_type: [], profile: [] };
      }
      f.density = f.cells.x.length
        ? CV.contours(f.cells.x, f.cells.y, b.meta.extent.cell)
        : { n: 0, levels: [], rings: [], reason: 'no cells' };
      delete f.cell_slice;
    }
    return {
      lineage: i, label: l.label, patient: l.patient,
      patient_id: b.names.patients[l.patient],
      receptor: l.receptor, n_cells: l.n_cells, frames: frames,
      axis: axis, fixed: fixed,
      min_cells_for_contour: CV.MIN_FOR_CONTOUR,
    };
  }

  function cset(b, q) {
    const layout = q.get('layout') || 'default';
    const uxy = uxyOf(b, layout);
    // Declared before the early return below, which reports them: a `const` referenced
    // above its declaration is a ReferenceError, not an undefined.
    const axis = q.get('axis') || 'time';
    let fixed = parseInt(q.get('fixed') || '-1', 10);
    const pats = q.getAll('patient').map((v) => parseInt(v, 10));
    const lins = new Set(q.getAll('lineage').map((v) => parseInt(v, 10)));
    if (pats.length) {
      const ps = new Set(pats);
      const L = b.names.lineages;
      for (let i = 0; i < L.length; i++) if (ps.has(L[i].patient)) lins.add(i);
    }
    if (!lins.size) {
      return { patients: pats, lineages: [], frames: [], n_lineages: 0, axis: axis, fixed: fixed,
               note: 'no clones selected' };
    }

    const inUnit = new Uint8Array(b.n_units);
    let nUnitsTotal = 0;
    for (let k = 0; k < b.n_units; k++) {
      if (lins.has(b.u_lineage[k])) { inUnit[k] = 1; nUnitsTotal++; }
    }
    const inCell = new Uint8Array(b.n_cells);
    let nCellsTotal = 0;
    for (let c = 0; c < b.n_cells; c++) {
      if (lins.has(b.c_lineage[c])) { inCell[c] = 1; nCellsTotal++; }
    }

    // Where each lineage sits at each visit, so an arrow can span a visit the lineage
    // missed rather than being dropped at the gap. server.py builds the same map.
    const pos = new Map();
    for (let k = 0; k < b.n_units; k++) {
      if (!inUnit[k]) continue;
      const li = b.u_lineage[k];
      if (!pos.has(li)) pos.set(li, new Map());
      pos.get(li).set(b.u_visit[k], k);
    }

    /* Same two axes as the server. On the tissue axis the frames are sites at ONE
     * visit, the clone map holds still, and arrows are suppressed because an arrow
     * between two anatomic positions would assert a movement Clonotrace never modelled. */
    if (axis === 'tissue' && fixed < 0) fixed = busiestVisit(b, (c) => inCell[c]);
    const tissues = tissuesOf(b);
    const steps = [];
    if (axis === 'tissue') {
      for (let ti = 0; ti < tissues.length; ti++) steps.push([fixed, ti]);
    } else {
      const ti = (fixed >= 0 && fixed < tissues.length) ? fixed : -1;
      for (let vi = 0; vi < b.visits.length; vi++) steps.push([vi, ti]);
    }

    const frames = [];
    for (const step of steps) {
      const vi = step[0], ti = step[1];
      const uidx = [], cidx = [];
      for (let k = 0; k < b.n_units; k++) if (inUnit[k] && b.u_visit[k] === vi) uidx.push(k);
      for (let c = 0; c < b.n_cells; c++) {
        if (!inCell[c] || b.c_visit[c] !== vi) continue;
        if (ti >= 0 && b.c_tissue[c] !== ti) continue;
        cidx.push(c);
      }
      if (axis === 'tissue' && cidx.length === 0) continue;

      // Contours from EVERY cell, before any subsampling, so the shape being read is
      // never the shape of a sample.
      let dmap;
      if (cidx.length) {
        const dx = new Float64Array(cidx.length), dy = new Float64Array(cidx.length);
        for (let t = 0; t < cidx.length; t++) {
          dx[t] = b.cxy[2 * cidx[t]]; dy[t] = b.cxy[2 * cidx[t] + 1];
        }
        dmap = CV.contours(dx, dy, b.meta.extent.cell);
      } else {
        dmap = { n: 0, levels: [], rings: [], reason: 'no cells' };
      }

      let shown = cidx, sub = false;
      if (cidx.length > MAX_CELLS_PER_FRAME) {
        // rng.choice without replacement, then sorted -- the same shape of operation as
        // server.py. The particular sample differs (see the note on mulberry32 above);
        // what matters is that it is a fixed one, so stepping back to a frame draws the
        // same cells rather than a new sample each time.
        const pool = shuffleInPlace(cidx.slice(), 12345);
        shown = pool.slice(0, MAX_CELLS_PER_FRAME).sort((p, r) => p - r);
        sub = true;
      }

      const arrows = [];
      if (axis !== 'tissue') pos.forEach(function (byv, li) {
        if (!byv.has(vi)) return;
        let pv = -1;
        byv.forEach(function (_, v) { if (v < vi && v > pv) pv = v; });
        if (pv < 0) return;
        const a = byv.get(pv), z = byv.get(vi);
        arrows.push({
          lineage: li,
          x1: uxy[2 * a], y1: uxy[2 * a + 1],
          x2: uxy[2 * z], y2: uxy[2 * z + 1],
          from_visit: pv, p1: b.u_profile[a], p2: b.u_profile[z],
          skipped: vi - pv > 1,
        });
      });

      const uniqLin = new Set();
      for (const k of uidx) uniqLin.add(b.u_lineage[k]);
      const sx = new Array(shown.length), sy = new Array(shown.length);
      for (let t = 0; t < shown.length; t++) {
        sx[t] = b.cxy[2 * shown[t]]; sy[t] = b.cxy[2 * shown[t] + 1];
      }

      const fr = {
        visit: vi,
        code: axis === 'tissue' ? tissues[ti] : b.visits[vi],
        label: axis === 'tissue' ? tissues[ti] : b.meta.visits[vi].label,
        n_units: uidx.length, n_cells: cidx.length,
        n_lineages: uniqLin.size,
        units: {
          x: uidx.map((k) => uxy[2 * k]),
          y: uidx.map((k) => uxy[2 * k + 1]),
          profile: uidx.map((k) => b.u_profile[k]),
          lineage: uidx.map((k) => b.u_lineage[k]),
          n_cells: uidx.map((k) => b.u_n[k]),
          ptrc: uidx.map((k) => b.u_ptrc[k]),
        },
        cells: { x: sx, y: sy },
        cells_shown: shown.length, subsampled: sub,
        density: dmap,
        arrows: arrows,
      };
      // Key order matters only to a diff, but PRESENCE matters to the parity suite and
      // to any consumer: server.py attaches these whenever a tissue is in play, so the
      // offline copy must too.
      if (ti >= 0) { fr.tissue = ti; fr.tissue_label = tissues[ti]; }
      frames.push(fr);
    }

    return {
      patients: pats,
      patient_ids: pats.map((p) => b.names.patients[p]),
      lineages_selected: q.getAll('lineage').map((v) => parseInt(v, 10)).sort((a, c) => a - c),
      axis: axis, fixed: fixed,
      n_lineages: lins.size,
      n_units: nUnitsTotal,
      n_cells: nCellsTotal,
      frames: frames,
      min_cells_for_contour: CV.MIN_FOR_CONTOUR,
      max_cells_per_frame: MAX_CELLS_PER_FRAME,
    };
  }

  /* ---- dispatch ------------------------------------------------------------ */

  const ROUTES = [
    ['/api/datasets', 'datasets', 0],
    ['/api/dataset/', 'dataset', 1],
    ['/api/centroids/', 'centroids', 1],
    ['/api/clones/', 'clones', 1],
    ['/api/set/', 'cset', 1],
    ['/api/lineages/', 'lineages', 1],
    ['/api/lineage/', 'lineage', 2],
  ];

  CV.route = function (url) {
    const qi = url.indexOf('?');
    const path = qi < 0 ? url : url.slice(0, qi);
    const q = new URLSearchParams(qi < 0 ? '' : url.slice(qi + 1));

    for (const [prefix, fn, nargs] of ROUTES) {
      if (nargs === 0 && path === prefix) return Promise.resolve(datasets());
      if (nargs && path.indexOf(prefix) === 0) {
        const rest = path.slice(prefix.length).split('/').filter((s) => s.length);
        if (rest.length !== nargs) continue;
        const name = rest[0];
        return CV.loadDataset(name).then(function (b) {
          switch (fn) {
            case 'dataset': return dataset(b);
            case 'centroids': return centroids(b, q);
            case 'clones': return clones(b, q);
            case 'cset': return cset(b, q);
            case 'lineages': return lineages(b, q);
            case 'lineage': return lineage(b, rest[1], q);
            default: throw new Error('unrouted: ' + fn);
          }
        });
      }
    }
    return Promise.reject(new Error('no such route: ' + path));
  };

  CV._api = { median, centroidsOf, lineageFrames, shuffleInPlace, mulberry32 };
}(window.CV));
